# NPS-ActiveSpace — Code conventions for humans and LLM agents

NPS-ActiveSpace is a scientific Python package for modeling the geographic extent of audibility from a given point using sound pressure level and causal data. It also includes related library functionality and visualization capabilities.

See [README.md](README.md) for detailed background and installation instructions. 

See [CONTRIBUTING.md](CONTRIBUTING.md) for PR workflow.

## Repo layout

- **Source**: `nps_active_space/` at the repo root. The codebase includes general library utils as well as workflow scripts (under `scripts/`).
- **Tests**: `tests/` — mirrors source paths; test file names mirror source modules
- **Config**: ConfigParser `.config` files in `nps_active_space/config/`; copy `template.config` to `<environment>.config` and pass `-e <environment>` to scripts.
- **Scripts**: run from the repo root, e.g. `python nps_active_space/scripts/run_ground_truthing.py -e DENA ...`. See [nps_active_space/scripts/README.md](nps_active_space/scripts/README.md) for more details on each script + example invocations.

## General Code Conventions
- Scientific logic lives in library modules (`active_space/`, `ground_truthing/`, `utils/`, `validation/`). Scripts handle argparse, config, then pass data and paths into library functions.
- Prefer concise, single-responsibility functions, classes, and modules. If one piece of logic gets too large, break it up into multiple clearly-named smaller ones.
- Don't leave comments that just restate the code. They should be used for scientific concepts, non-obvious geospatial operations, or to explain decisions.
- Units should be clearly specified in variable names (e.g., `altitude_m`, `EVENT_GAP_SECONDS`, `elapsed_seconds`). Unit conversions should live in clearly named functions (e.g. `feet_to_meters`, `seconds_to_hours`)


## Python Code Conventions
- Target Python 3.12.
- Use modern type hints for all function signatures: `list[X]`, `dict[str, X]`, `X | None`. Not `typing.List`, `typing.Optional`.
- Prefer `match`/`case` for handling enum values or exhaustive case matching.
- Name DataFrames for what they hold (e.g. `sorted_mxak_points`, `flight_tracks`), not just `df`.
- Use `pathlib` and ensure cross-platform path support for any code dealing with file system paths.
- Timestamps are currently stored as tz-naive. Data parsers should return UTC-naive. See our [UTC standardization issue](https://github.com/dbetchkal/NPS-ActiveSpace/issues/96) for desired future state.

## Testing

- `pytest` from the repo root with an activated venv.
- Group related tests into classes (`class TestX:`).
- DataFrames: use pandas or geopandas `testing.assert_frame_equal` against an expected result for test assertions on dataframe results. For narrow tests, assert a single field or created an expected dataframe with a subset of columns.
- New logic should ship with tests. If something is difficult to test in a unit test (physics model subprocess, GUI, external binaries), mention it in the PR description explaining why. 
- Document any manual testing/verification that you do in PR descriptions.

## Git
- Keep commits logically scoped.
- After a squash-merge, rebase dependent PRs with the `--onto` flag to minimize merge conflicts: `git rebase --onto <new-base> <old-base> <branch>`.

## Additional agent instructions

- Activate `.venv` before running Python commands.
- Run `pytest` after code changes. Don't mark work complete without passing tests.
- Keep changes small and scoped. Don't refactor unrelated code or rename/move files without maintainer sign-off.
