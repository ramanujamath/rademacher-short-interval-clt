"""Run every verification script and report a combined status.

Usage:  python3 run_all.py [--quick]
--quick lowers the Monte Carlo sample size (1000 instead of 4000).
"""
import subprocess
import sys
import time
from pathlib import Path

here = Path(__file__).resolve().parent
quick = "--quick" in sys.argv
scripts = [
    ["check_kernel.py"],
    ["check_local_moments.py"],
    ["check_prime_sums.py"],
    ["check_combinatorics.py"],
    ["check_transfer_identities.py"],
    ["check_asymptotic_bookkeeping.py"],
    ["monte_carlo_clt.py", "1000" if quick else "4000"],
]
failed = []
for cmd in scripts:
    t0 = time.time()
    r = subprocess.run([sys.executable, str(here / cmd[0]), *cmd[1:]], cwd=here)
    print(f"--> {cmd[0]}: {'OK' if r.returncode == 0 else 'FAILED'} ({time.time() - t0:.0f}s)\n")
    if r.returncode:
        failed.append(cmd[0])
print("ALL CHECKS PASSED" if not failed else f"FAILED: {', '.join(failed)}")
sys.exit(1 if failed else 0)
