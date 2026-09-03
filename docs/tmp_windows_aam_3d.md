# Windows: 3D active space AAM vs NMSim (DENATRLA 2025)

Native Windows run of the **full 3D product path** on `example_data` DENATRLA:

```text
generate_3d_active_space.py → generate_active_space_batch.py → fit_3d_active_space.py
```

Same pipeline NMSim production uses; AAM passes `--model aam` on each layer command.

## Setup

Same as [tmp_windows_aam_1500m.md](tmp_windows_aam_1500m.md): `AK.config` from `AK_example.config`, set `[project] nmsim` and `[project] aam`.

## Run

```bat
docs\tmp_windows_aam_3d.bat
docs\tmp_windows_aam_3d.bat aam
docs\tmp_windows_aam_3d.bat nmsim
```

Default **both** models. Use `aam` only if NMSim 3D layers already exist locally (common on dev Macs — not committed to git).

### Quick example smoke (recommended first on Windows)

**~3 layers**, macparity mesh (density 10, omni `+000`), committed annotations:

```bat
docs\tmp_windows_aam_3d_example.bat
docs\tmp_windows_aam_3d_example.bat aam
```

Layers **1200, 1500, 1800 m**. Metadata: [tmp_platform_aam_3d_example_job.json](tmp_platform_aam_3d_example_job.json). Mac Docker: `docs/tmp_docker_aam_3d_example.sh`.

### Full standard run

| Setting | Value |
|---------|--------|
| Altitudes | **900, 1200, 1500, 1800, 2100 m** (300 m step) |
| Density | **48** (`DEFAULT_SRC_PT_DENSITY`; same as `tmp_windows_aam_1500m.bat` full) |
| Omni | **0–2 dB** (0, +0.5, +1, +1.5, +2) |
| Headings | 0 |
| Annotations | `example_data\site_projects\DENATRLA\DENATRLA2025_saved_annotations.geojson` (bat passes basename; generate joins site dir) |

### What the full job uses

Metadata: [tmp_platform_aam_3d_job.json](tmp_platform_aam_3d_job.json).

Mac Docker equivalent: `docs/tmp_docker_aam_3d.sh`.

### Runtime (Windows host)

Mac-parity single layer (`docs\tmp_windows_aam_macparity.bat`: 1500 m, density 10, omni `+000`) was about **204 s AAM** and **108 s NMSim** on the Windows host.

This 3D job is **five layers**, **~5× omni gains**, and a **much denser mesh** — expect **hours**, not minutes. For a quick pipeline check only, temporarily edit `EXTRA=` in the bat to macparity settings (`--density 10 --omni-max 0`).

## Outputs

Under `example_data\site_projects\DENATRLA\`:

| Path | Role |
|------|------|
| `Output_Data\{aam\|nmsim}\ACTIVESPACES\DENATRLA2025_{alt}m\` | Per-layer GeoJSON (`O_+000`, …) |
| `Output_Data\{model}\PRECISION_RECALL\PrecisionRecallPlot_3d_DENATRLA2025_f1p0.png` | 3D gain fit PR plot |
| `Output_Data\AMBIENCE\DENATRLA2025_ambience.pkl` | Shared NVSPL cache (first run) |
| `DENATRLA2025_commands.txt` | Batch command file (overwritten each model run) |

Project **`example_data\site_projects\fits.csv`** gets one row per model (`Model` = `aam` / `nmsim`).

## Viz

After `fit_3d_active_space.py`, viz can use the fitted gain (omit `-g`). Basename is OK for `--annotation-file` (resolved under the site dir); full path also works:

```bat
python nps_active_space\scripts\viz.py DENATRLA2025 -e AK -s -a --compare --terraced ^
  --annotation-file DENATRLA2025_saved_annotations.geojson
```

## Existing NMSim 3D layers

The repo does **not** ship generated `ACTIVESPACES` trees (gitignored). Check your site dir:

```bat
dir example_data\site_projects\DENATRLA\Output_Data\nmsim\ACTIVESPACES
dir example_data\site_projects\DENATRLA\Output_Data\aam\ACTIVESPACES
```

If NMSim already has `DENATRLA2025_900m` … `_2100m`, run `docs\tmp_windows_aam_3d.bat aam` to build only the AAM side, then compare in viz.
