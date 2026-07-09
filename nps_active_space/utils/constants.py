"""Shared constants for the NPS-ActiveSpace package."""
import os

# os.name on Windows — see https://docs.python.org/3/library/os.html#os.name
OS_NAME_WINDOWS = "nt"
IS_WINDOWS = os.name == OS_NAME_WINDOWS

# TODO: put global unit conversion constants here (e.g. FEET_TO_METERS from project_setup.py)
