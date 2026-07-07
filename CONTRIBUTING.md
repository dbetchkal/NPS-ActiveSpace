# Contributing to NPS-ActiveSpace

Thank you for your interest in contributing to this project. Whether you're fixing
a bug, adding a feature, or improving documentation, your help is welcome.

## Getting started

1. Fork and clone the repository:

   ```
   git clone https://github.com/<your-username>/NPS-ActiveSpace.git
   cd NPS-ActiveSpace
   ```

2. Create a virtual environment (Python 3.12+):

   ```
   python -m venv .venv
   source .venv/bin/activate
   ```

3. Install system dependencies. On Ubuntu/Debian:

   ```
   sudo apt-get install gdal-bin libgdal-dev
   ```

   On macOS with Homebrew: `brew install gdal`. The GDAL Python binding version
   must match your system GDAL version.

4. Install the project and its dependencies:

   ```
   pip install --upgrade pip
   pip install -r requirements.txt
   pip install -e .
   pip install pytest
   ```

5. Run the test suite to make sure everything works:

   ```
   pytest tests/ -v
   ```

## Making changes

- Create a feature branch from `main` (e.g. `feature/my-improvement` or `fix/issue-description`).
- Write tests for new functionality when possible. Tests live in the `tests/` directory.
- Keep commits focused. One logical change per commit is easier to review.

## Pull requests

When you're ready to submit your work:

- Give the PR a clear, descriptive title.
- Explain what the change does and why it's needed in the PR description.
- Make sure CI passes (the test suite runs automatically on PRs to `main`).
- If your change relates to an open issue, reference it in the description.

Smaller, well-scoped PRs are easier to review and more likely to be merged quickly.

## Code style

- Follow the patterns you see in the existing codebase.
- Python 3.12+ features are fine to use.
- Add type hints where practical, especially for function signatures.
- Configuration files use the ConfigParser `.config` format (see `nps_active_space/config/`
  for examples).

## License

This project is in the public domain within the United States under
[CC0 1.0 Universal](LICENSE.md). All contributions will be released under the
same terms. By submitting a pull request, you agree to this waiver of copyright
interest.
