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
import os
import queue
import subprocess
import threading
import time

DEFAULT_FREE_MB = 4000          # a card with less than this in use counts as free
POLL_SECONDS = 15


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


class GPUPool:
    """Hand out device ids, waiting for one that is genuinely idle."""

    def __init__(self, gpus, free_mb=DEFAULT_FREE_MB, check_external=True):
        self.free_mb = free_mb
        self.check_external = check_external
        self._q = queue.Queue()
        for g in gpus:
            self._q.put(int(g))
        self._lock = threading.Lock()

    def acquire(self, timeout=None):
        """Block until a device is free in the pool AND quiet on the card."""
        deadline = None if timeout is None else time.time() + timeout
        while True:
            g = self._q.get() if deadline is None else self._q.get(
                timeout=max(0.1, deadline - time.time()))
            if not self.check_external:
                return g
            used = gpu_memory_used()
            if not used or used.get(g, 0) < self.free_mb:
                return g
            # someone else is on this card: put it back and try another
            self._q.put(g)
            time.sleep(POLL_SECONDS)

    def release(self, g):
        self._q.put(int(g))

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
