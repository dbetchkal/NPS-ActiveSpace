import sys
from pathlib import Path

_setup_tests = Path(__file__).resolve().parent
if str(_setup_tests) not in sys.path:
    sys.path.insert(0, str(_setup_tests))
