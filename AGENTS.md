# NPS-ActiveSpace — Code Conventions

NPS-ActiveSpace is a scientific Python package for modeling geographic extent of audibility from a given point using sound pressure level and causal data. It also includes related library functionality and visualization capabilities.

See [README.md](README.md) for detailed background and installation instructions. 

See [CONTRIBUTING.md](CONTRIBUTING.md) for PR workflow.

> Activate `.venv` before running Python commands.

## Repo layout

- **Source**: `nps_active_space/` at the repo root. The codebase includes general library utils as well as workflow scripts (under `scripts/`).
- **Tests**: `tests/` — mirrors source paths; test file names mirror source modules
- **Config**: ConfigParser `.config` files in `nps_active_space/config/`; copy `template.config` to `<environment>.config` and pass `-e <environment>` to scripts.
- **Scripts**: run from the repo root, e.g. `python nps_active_space/scripts/run_ground_truthing.py -e DENA ...`. See [nps_active_space/scripts/README.md](nps_active_space/scripts/README.md) for more details on each script + example invocations.

## General Code Conventions
- Scientific logic lives in library modules (`active_space/`, `ground_truthing/`, `utils/`, `validation/`). Scripts handle argparse, config, and I/O, then pass data and paths into library functions.
- Prefer concise, single-responsibility functions. If one function or class gets too large, break it into multiple smaller ones.
- Use descriptive variable and function names.
- Comments explain the unclear *why* and should be used for scientific concepts, potentially confusing geospatial operations, or non-obvious trade-offs. Don't just restate the code.
- Units should be clearly specified in variable names (e.g., `altitude_m`, `EVENT_GAP_SECONDS`, `elapsed_seconds`)


## Python Code Conventions
- Target Python 3.12. 
- Use type hints for all new code and prefer modern type hints only: `list[X]`, `dict[str, X]`, `X | None`. Not `typing.List`, `typing.Optional`.
- Prefer `match`/`case` in named library functions for handling enums or exhaustive case matching.
- Name DataFrames for what they hold (e.g. `mxak_points`, `flight_tracks`), not `df`.

## Testing

- `pytest` from the repo root with an activated venv.
- Group related tests into classes (`class TestX:`).
- DataFrames: use pandas or geopandas `testing.assert_frame_equal` for test assertions on dataframe results. Dataclasses/Pydantic: assert equality against an expected instance. For narrow tests, assert a single field or geometry check.
- Algorithm logic: synthetic/minimal unit tests. Real-format I/O: parametrized regression against bundled fixtures.
- New logic ships with tests. If something is difficult to test in a unit test (physics model subprocess, GUI, external binaries), mention it in the PR description and document any manual testing/verification that you do.

## Git

- Propose file renames, subdir moves, or structural refactors for maintainer sign-off before implementing.
- Keep commits logically scoped.
- After a squash-merge, rebase dependent PRs: `git rebase --onto <new-base> <old-base> <branch>`.
