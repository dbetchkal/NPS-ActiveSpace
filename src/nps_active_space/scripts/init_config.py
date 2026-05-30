"""Create a user config file from the bundled template."""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

import nps_active_space.config as cfg


def main() -> None:
    parser = ArgumentParser(
        description="Create an environment config TOML from the bundled template."
    )
    parser.add_argument(
        "-e",
        "--environment",
        required=True,
        help="Config name (writes {environment}.toml, e.g. DENA).",
    )
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=None,
        help=(
            "Destination folder (default: NPS_ACTIVE_SPACE_CONFIG_DIR or "
            "~/.nps_active_space)."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing config file.",
    )
    args = parser.parse_args()

    dest = cfg.init_config(
        args.environment,
        config_dir=args.config_dir,
        force=args.force,
    )
    print(f"Wrote config: {dest}")
    print("Edit paths and credentials, then run pipeline scripts with -e", args.environment)
    if hint := cfg.config_setup_hint(dest.parent):
        print(hint)


if __name__ == "__main__":
    main()
