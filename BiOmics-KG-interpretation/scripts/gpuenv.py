"""Re-exec the unified brick CLI with the cluster's CUDA 11.8 runtime."""

from __future__ import annotations

import os
import sys
from pathlib import Path


DEFAULT_LIBRARY_DIRS = (
    "/home/HPCBase/compilers/gcc/9.3.0/lib64",
    "/home/HPCBase/compilers/cuda/11.8.0/lib64",
    "/home/HPCBase/compilers/cuda/11.8.0/extras/CUPTI/lib64",
    "/home/HPCBase/libs/cudnn/8.4.0.27_cuda11.x/lib",
    "/home/HPCBase/libs/openblas/0.3.25_gcc9.3.0/lib64",
)


def _bundled_sklearn_openmp():
    candidates = sorted(
        (Path(sys.prefix) / "lib" / f"python{sys.version_info.major}.{sys.version_info.minor}" / "site-packages" / "scikit_learn.libs").glob("libgomp*.so*")
    )
    return str(candidates[0]) if candidates else ""


def _runtime_env():
    env = os.environ.copy()
    configured = env.get("BRICK_GPU_LIBRARY_PATH")
    library_dirs = configured.split(os.pathsep) if configured else list(DEFAULT_LIBRARY_DIRS)
    missing = [path for path in library_dirs if not Path(path).is_dir()]
    if missing:
        raise RuntimeError("Missing GPU runtime directories: " + ", ".join(missing))
    existing = [value for value in env.get("LD_LIBRARY_PATH", "").split(os.pathsep) if value]
    env["LD_LIBRARY_PATH"] = os.pathsep.join([*library_dirs, *existing])
    preload = env.get("BRICK_GPU_LD_PRELOAD") or _bundled_sklearn_openmp()
    if not preload or not Path(preload).exists():
        raise RuntimeError(
            "Could not locate the scikit-learn OpenMP runtime required on this aarch64 cluster; "
            "set BRICK_GPU_LD_PRELOAD explicitly"
        )
    existing_preload = [value for value in env.get("LD_PRELOAD", "").split() if value]
    env["LD_PRELOAD"] = " ".join([preload, *existing_preload])
    env["BRICK_GPU_ENV_READY"] = "1"
    return env


def main():
    if os.environ.get("BRICK_GPU_ENV_READY") != "1":
        env = _runtime_env()
        code = "from cli import main; main()"
        os.execve(sys.executable, [sys.executable, "-c", code, *sys.argv[1:]], env)
    from cli import main as brick_main

    brick_main()


if __name__ == "__main__":
    main()
