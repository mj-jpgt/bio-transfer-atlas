import subprocess
import sys
from pathlib import Path


def run(script: Path):
    print(f"\n[gate] {script}")
    r = subprocess.run([sys.executable, str(script)])
    if r.returncode != 0:
        raise SystemExit(r.returncode)


if __name__ == "__main__":
    gates = [
        Path("tests/gate_finngen.py"),
        Path("tests/gate_bbj.py"),
        Path("tests/gate_multisource_labels.py"),
        Path("tests/gate_constraint.py"),
        Path("tests/gate_eval.py"),
    ]
    for g in gates:
        if g.exists():
            run(g)
    print("\nAll gate scripts passed.")
