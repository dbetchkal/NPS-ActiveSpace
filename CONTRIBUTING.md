# Contributing

## Setup

Install as described in the [README](README.md#installation). Use `pip install -e ".[dev]"` to include test tools.

## Workflow

1. Fork and create a branch — `feature/description` or `fix/description`.
2. New logic should ship with tests. See [Style Guide](AGENTS.md) for more details.
3. Ensure existing and new tests pass before marking a PR ready for review. GitHub actions will run against PRs pointing at main.
4. Open a pull request against `main` with a clear title and description of what changed and why. Reference any related issues.

Prefer small, focused PRs with a concise, clear description.

## Code style

See [`AGENTS.md`](AGENTS.md) for Python conventions, naming, type hints, testing patterns, and git workflow.

## License

This project is [public domain (CC0)](LICENSE.md). By submitting a pull request you waive any copyright interest in your contribution.
