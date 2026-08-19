"""Validate an environment config file."""

import argparse
import sys

import nps_active_space.utils.config as cfg

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Validate a NPS-ActiveSpace environment config file.",
    )
    parser.add_argument(
        "-e",
        "--environment",
        required=True,
        help="Configuration environment to validate (e.g. DENA_example).",
    )
    args = parser.parse_args()

    try:
        cfg.initialize(args.environment, validate_config=True)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)

    print(f"Configuration {args.environment!r} is valid.")
