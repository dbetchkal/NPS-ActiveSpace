#!/usr/bin/env bash
# Install the GDAL Python binding for macOS/Linux (source build against system libgdal),
# then optionally install this repo in editable mode. Windows uses a wheel in pyproject.toml.
#
# Prerequisite: system GDAL (e.g. brew install gdal; apt install libgdal-dev python3.12-dev).
# Run from the repository root with an activated venv.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Keep in sync with pyproject.toml [project.dependencies] numpy pin.
NUMPY_VERSION="2.3.1"

EXTRAS="dev"
SKIP_EDITABLE=false
NO_VENV_CHECK=false
SET_INCLUDE_PATHS=false

usage() {
  cat <<'EOF'
Usage: install_scripts/install_gdal_python_binding.sh [OPTIONS]

Build and pip-install GDAL==$(gdal-config --version) (PyPI has no Linux/macOS wheels),
then pip install -e ".[<extras>]" from the repo root.

Options:
  --extras EXTRAS       Optional dependency group(s), comma-separated (default: dev)
  --gdal-only           Skip editable install of NPS-ActiveSpace
  --no-venv-check       Allow running without VIRTUAL_ENV (e.g. GitHub Actions)
  --set-include-paths   Export C/C++ include paths for Debian-style /usr/include/gdal
  -h, --help            Show this help

Environment:
  GDAL_CONFIG           Path to gdal-config if not on PATH
  GDAL_VERSION          Override version (default: gdal-config --version)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --extras)
      EXTRAS="${2:?--extras requires a value}"
      shift 2
      ;;
    --gdal-only)
      SKIP_EDITABLE=true
      shift
      ;;
    --no-venv-check)
      NO_VENV_CHECK=true
      shift
      ;;
    --set-include-paths)
      SET_INCLUDE_PATHS=true
      shift
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ "${NO_VENV_CHECK}" != true && -z "${VIRTUAL_ENV:-}" ]]; then
  cat >&2 <<'EOF'
No active virtual environment (VIRTUAL_ENV is unset).
Create and activate a venv, then re-run this script:
  python3.12 -m venv .venv && source .venv/bin/activate
Or pass --no-venv-check when installing into a dedicated Python (e.g. CI).
EOF
  exit 1
fi

GDAL_CONFIG="${GDAL_CONFIG:-gdal-config}"
if ! command -v "${GDAL_CONFIG}" >/dev/null 2>&1; then
  cat >&2 <<EOF
${GDAL_CONFIG} not found. Install system GDAL first (README Installation section).
  macOS: brew install gdal  (ensure \$(brew --prefix gdal)/bin is on PATH)
  Debian/Ubuntu: sudo apt-get install gdal-bin libgdal-dev python3.12-dev build-essential
  Or set GDAL_CONFIG to the full path of gdal-config.
EOF
  exit 1
fi

GDAL_VERSION="${GDAL_VERSION:-$("${GDAL_CONFIG}" --version)}"
echo "Using GDAL Python binding version ${GDAL_VERSION} (from ${GDAL_CONFIG})"

cd "${REPO_ROOT}"

python -m pip install --upgrade pip wheel
python -m pip install "setuptools<84"
python -m pip install "numpy==${NUMPY_VERSION}"

install_gdal() {
  if [[ "${SET_INCLUDE_PATHS}" == true ]]; then
    export CPLUS_INCLUDE_PATH="/usr/include/gdal${CPLUS_INCLUDE_PATH:+:${CPLUS_INCLUDE_PATH}}"
    export C_INCLUDE_PATH="/usr/include/gdal${C_INCLUDE_PATH:+:${C_INCLUDE_PATH}}"
  fi
  python -m pip install "GDAL==${GDAL_VERSION}" --no-build-isolation
}

install_gdal

if [[ "${SKIP_EDITABLE}" != true ]]; then
  python -m pip install -e ".[${EXTRAS}]"
fi

echo "GDAL Python binding install complete."
