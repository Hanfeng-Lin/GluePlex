"""Launch the pass-2 Boltz run scripts emitted by 4.1.prepare_second_pass.py.

Reads pass2_manifest.txt and runs each wrapper python script in parallel.
Parallelism is bounded by the number of unique GPUs referenced in the scripts
(each run sets its own CUDA_VISIBLE_DEVICES), matching the 3.1 pattern.
"""

import argparse
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _gpu_pool import GPUPool, slot


def extract_gpu(script_path: Path) -> int:
    m = re.search(r'CUDA_VISIBLE_DEVICES"\]\s*=\s*str\((\d+)\)', script_path.read_text())
    return int(m.group(1)) if m else 0


def run_one(script_path: Path, cwd: Path, pool: GPUPool) -> int:
    """Run one pass-2 job on whichever card is free.

    The GPU baked into the script by 4.1 is only a fallback now: the scripts use
    os.environ.setdefault, so the device pinned here by the launcher wins. Fixed
    assignment left cards idle whenever another pipeline was sharing the host.
    """
    with slot(pool) as gpu:
        print(f"[START] {script_path.name}  (GPU {gpu})")
        try:
            pool.run([sys.executable, str(script_path)], gpu,
                     check=True, cwd=str(cwd))
            print(f"[ OK ] {script_path.name}")
            return 0
        except subprocess.CalledProcessError as e:
            print(f"[FAIL] {script_path.name}: {e}", file=sys.stderr)
            return e.returncode or 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step4_dir", required=True,
                        help="step4_boltz directory containing pass2_manifest.txt and run scripts.")
    args = parser.parse_args()

    step4_dir = Path(args.step4_dir).resolve()
    manifest = step4_dir / "pass2_manifest.txt"
    if not manifest.exists():
        print(f"Error: manifest {manifest} not found. Run 4.1 first.", file=sys.stderr)
        sys.exit(1)

    scripts = [Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip()]
    if not scripts:
        print("Manifest is empty — nothing to run.", file=sys.stderr)
        sys.exit(1)

    # The GPUs the scripts were generated against define the pool; which job
    # lands on which card is decided at run time.
    unique_gpus = sorted({extract_gpu(s) for s in scripts})
    gpu_pool = GPUPool(unique_gpus)
    max_workers = max(1, len(unique_gpus))
    print(f"Launching {len(scripts)} jobs across {max_workers} worker(s); "
          f"GPU pool: {unique_gpus}")

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        results = list(ex.map(lambda s: run_one(s, step4_dir, gpu_pool), scripts))
    failures = sum(1 for r in results if r != 0)

    if failures:
        print(f"{failures}/{len(scripts)} pass-2 jobs failed.", file=sys.stderr)
        sys.exit(1)
    print(f"All {len(scripts)} pass-2 jobs completed.")


if __name__ == "__main__":
    main()
