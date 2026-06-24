import sys
from pathlib import Path

# Allow `from helpers import ...` in tests/ground_truthing/
_gt_tests = Path(__file__).resolve().parent / "ground_truthing"
if str(_gt_tests) not in sys.path:
    sys.path.insert(0, str(_gt_tests))
