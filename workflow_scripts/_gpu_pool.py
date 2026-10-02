#!/usr/bin/env python
"""A GPU pool that hands out whichever device is actually free.

Replaces static round-robin (`gpus[i % len(gpus)]`), which had two problems:

1. **It raced on the environment.** Each worker did
   `os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)` before launching its
   subprocess, but threads share one environ, so whichever thread wrote last
   won and several jobs landed on the same card while others sat idle. The
   device must be passed to the child through `subprocess(env=...)` instead.

2. **It could not see other jobs.** Round-robin assumes it owns the machine.
   When two pipelines run side by side (say a truncation arm and a docking arm)
   both hand out device 0 first and collide.

`acquire()` blocks until a device is free both in this pool and on the card
itself (checked with nvidia-smi), so concurrent pipelines interleave correctly.
"""
import errno
import fcntl
import os
import queue
import subprocess
import threading
import time

DEFAULT_FREE_MB = 4000          # a card with less than this in use counts as free
POLL_SECONDS = 15
LOCK_DIR = os.environ.get("GLUEPLEX_GPU_LOCK_DIR", "/tmp/glueplex_gpu_locks")


class _GpuLock:
    """Exclusive, cross-process claim on one device, via flock.

    An in-process queue plus an nvidia-smi check is not enough: two pipelines
    started in the same second each pop device 0 from their own queue, both see
    it idle, and both take it -- the check and the claim are not atomic. The
    lock makes them atomic, and it is released automatically if a process dies
    because flock is tied to the file descriptor.
    """

    def __init__(self, gpu):
        self.gpu = int(gpu)
        os.makedirs(LOCK_DIR, exist_ok=True)
        self.path = os.path.join(LOCK_DIR, f"gpu{self.gpu}.lock")
        self.fh = None

    def acquire(self):
        fh = open(self.path, "w")
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as e:
            fh.close()
            if e.errno in (errno.EACCES, errno.EAGAIN):
                return False
            raise
        fh.write(f"{os.getpid()}\n")
        fh.flush()
        self.fh = fh
        return True

    def release(self):
        if self.fh is not None:
            try:
                fcntl.flock(self.fh, fcntl.LOCK_UN)
            finally:
                self.fh.close()
                self.fh = None


def gpu_memory_used(samples=3, gap=1.5):
    """{index: peak MiB in use over `samples` readings}.

    A single reading is useless here: boltz starts a fresh process per model, so
    a card in active use reports 16 GB one second and 0.7 GB the next. Taking the
    peak over a few seconds distinguishes "between models" from "idle".
    """
    peak = {}
    for k in range(samples):
        for i, m in _read_once().items():
            peak[i] = max(peak.get(i, 0), m)
        if k < samples - 1:
            time.sleep(gap)
    return peak


def _read_once():
    """{index: MiB in use}; empty dict if nvidia-smi is unavailable."""
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=index,memory.used", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=20)
        used = {}
        for line in out.stdout.strip().splitlines():
            i, m = (x.strip() for x in line.split(","))
            used[int(i)] = int(m)
        return used
    except Exception:
        return {}


def _gpus_from_proc(samples=3, gap=1.5):
    """{index} of GPUs a live process is pinned to, read from /proc.

    nvidia-smi on this host sometimes reports `memory.used` as 1 MiB for every
    card and an empty `--query-compute-apps` while jobs are plainly running
    (nvidia-dcgm wedged; exit code is still 0). That silently disables the
    external-process check below, because no card ever looks busy. /proc does
    not lie: a process pinned with CUDA_VISIBLE_DEVICES=<i> that also holds an
    open /dev/nvidia* handle is on that card.

    The nvidia-fd requirement matters -- wrapper processes (conda run, the
    pipeline driver) inherit CUDA_VISIBLE_DEVICES without touching a GPU, and
    counting those would back off from cards that are actually free.

    Sampled over a window and unioned for the same reason `gpu_memory_used`
    peaks over one: boltz starts a fresh process per model, so a single scan
    taken between two models sees an empty /proc and calls a busy card free.
    """
    busy = set()
    for k in range(samples):
        busy |= _scan_proc_once()
        if k < samples - 1:
            time.sleep(gap)
    return busy


def _busy_in_proc(gpu, samples=3, gap=1.5):
    """Is `gpu` in use according to /proc? Returns as soon as it is seen.

    Only the "looks free" answer needs the full window -- that is the one a
    gap between two boltz models can fake. Seeing the card in use once is
    already conclusive, so the common busy case costs a single scan.
    """
    for k in range(samples):
        if gpu in _scan_proc_once():
            return True
        if k < samples - 1:
            time.sleep(gap)
    return False


def _scan_proc_once():
    busy = set()
    try:
        pids = [d for d in os.listdir("/proc") if d.isdigit()]
    except OSError:
        return busy
    for pid in pids:
        try:
            with open(f"/proc/{pid}/environ", "rb") as f:
                env = f.read()
        except (OSError, PermissionError):
            continue
        idx = None
        for kv in env.split(b"\0"):
            if kv.startswith(b"CUDA_VISIBLE_DEVICES="):
                v = kv.split(b"=", 1)[1].decode("ascii", "ignore").strip()
                if v.isdigit():
                    idx = int(v)
                break
        if idx is None:
            continue
        try:
            fds = os.listdir(f"/proc/{pid}/fd")
        except (OSError, PermissionError):
            continue
        for fd in fds:
            try:
                if os.readlink(f"/proc/{pid}/fd/{fd}").startswith("/dev/nvidia"):
                    busy.add(idx)
                    break
            except OSError:
                continue
    return busy


class GPUPool:
    """Hand out device ids, waiting for one that is genuinely idle."""

    def __init__(self, gpus, free_mb=DEFAULT_FREE_MB, check_external=True):
        self.free_mb = free_mb
        self.check_external = check_external
        self._q = queue.Queue()
        for g in gpus:
            self._q.put(int(g))
        self._lock = threading.Lock()
        self._held = {}

    def acquire(self, timeout=None):
        """Block until a device is free in this pool and claimable machine-wide."""
        deadline = None if timeout is None else time.time() + timeout
        while True:
            # The deadline has to be checked here, not only on the queue get:
            # a card rejected as busy is put back, so the queue is never empty
            # and `get` returns immediately forever. Before the external check
            # could actually fire this loop never spun, which hid the bug.
            if deadline is not None and time.time() >= deadline:
                raise queue.Empty(
                    f"no GPU in {sorted(self._held) or 'pool'} became free "
                    f"within {timeout}s")
            g = self._q.get() if deadline is None else self._q.get(
                timeout=max(0.1, deadline - time.time()))
            lock = _GpuLock(g)
            if not lock.acquire():
                # held by another pipeline on this host
                self._q.put(g)
                time.sleep(POLL_SECONDS)
                continue
            if self.check_external:
                used = gpu_memory_used()
                busy = used.get(g, 0) >= self.free_mb if used else False
                # nvidia-smi can report every card as idle while jobs run, so a
                # clean reading is not proof the card is free -- confirm via /proc
                # before trusting it.
                if not busy and (not used or max(used.values()) < self.free_mb):
                    busy = _busy_in_proc(g)
                if busy:
                    # busy with something that is not using these locks
                    lock.release()
                    self._q.put(g)
                    time.sleep(POLL_SECONDS)
                    continue
            with self._lock:
                self._held[g] = lock
            return g

    def release(self, g):
        g = int(g)
        with self._lock:
            lock = self._held.pop(g, None)
        if lock is not None:
            lock.release()
        self._q.put(g)

    def run(self, cmd, gpu, **kw):
        """subprocess.run with CUDA_VISIBLE_DEVICES set for the CHILD only."""
        env = dict(os.environ)
        env["CUDA_VISIBLE_DEVICES"] = str(gpu)
        return subprocess.run(cmd, env=env, **kw)


class _Slot:
    def __init__(self, pool):
        self.pool = pool
        self.gpu = None

    def __enter__(self):
        self.gpu = self.pool.acquire()
        return self.gpu

    def __exit__(self, *exc):
        if self.gpu is not None:
            self.pool.release(self.gpu)
        return False


def slot(pool):
    """`with slot(pool) as gpu:` — acquire for the duration of the block."""
    return _Slot(pool)
