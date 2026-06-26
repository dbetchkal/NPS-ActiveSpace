# Example data

Small, real-format samples for offline unit tests, local dev (`GLBA_example.config`, `DENA_example.config`), and manual debugging.

Examples days:
* The Glacier Bay example config points at May 25th, 2024 data for AIS vessel tracks and sound pressure level (NVSPL) at the Lower South Tidal (GLBALSTL) acoustic monitoring site.
* The Denali example config points at June 23rd, 2025 data for ADS-B flight tracks and sound pressure level (NVSPL) at the Triple Lakes Trail (DENATRLA) acoustic monitoring site.

Paths in the example configs should parse on Mac/Linux or Windows. However, the paths are relative, so scripts need to be run from the repository root.

**Clock conventions:** MXAK AIS timestamps are UTC-naive; NVSPL and ship-visit times are site-local naive (e.g. GLBALSTL ≈ `America/Juneau`).

## Layout

| Path | Purpose |
|------|---------|
| `nvspl_archive/` | iyore NVSPL archive (`.structure.txt` + deployment folders) |
| `site_projects/` | Per-site workspaces (`GLBALSTL/`, `DENATRLA/`) — DEM, study area, `.sit` |
| `MXAK-AIS-GLBA/` | Daily MXAK AIS CSVs for Glacier Bay |
| `ADS-B/healy_repeater/` | ADS-B TSV samples (Healy Repeater in Denali) |
| `faa/` | FAA type corrections + trimmed releasable-aircraft DB sample for DENATRLA example day |
| `GLBALSTL_ship_visits.csv` | Extra file with closest passes of ships by the GLBALSTL site (for general reference and AIS/NVSPL alignment tests) |

### GLBALSTL 2024-05-24 slice (~50 MB total)

| Path | Size (approx) | Notes |
|------|---------------|-------|
| `nvspl_archive/2024 GLBALSTL Lower South Tidal Inlet/01 DATA/NVSPL/` | 23 MB | 24 hourly `NVSPL_GLBALSTL_2024_05_24_*.txt` |
| `site_projects/GLBALSTL/` | 3.3 MB | `elevation_m_nad83_utm8.tif`, study area, `GLBALSTL2024.sit` |
| `MXAK-AIS-GLBA/MXAK-AIS-GLBA-20240524.csv` | 2.5 MB | Full-day AIS for alignment with NVSPL |

### DENATRLA 2025-06-23 slice (~25 MB new; shared archive + ADS-B)

| Path | Size (approx) | Notes |
|------|---------------|-------|
| `nvspl_archive/2025 DENATRLA Triple Lakes/01 DATA/NVSPL/` | 23 MB | 24 hourly `NVSPL_DENATRLA_2025_06_23_*.txt` |
| `site_projects/DENATRLA/` | 1.7 MB | `elevation_m_nad83_utm6.tif`, study area, `DENATRLA2025.sit` |
| `ADS-B/healy_repeater/20250623.TSV` | 4.4 MB | Shared with `test_adsb.py` — same calendar day as NVSPL |
| `faa/MASTER_DENATRLA_20250623_sample.txt` | 6 KB | 27 FAA registry rows for ICAOs on 2025-06-23 in study area (2 ICAOs not in US registry) |

## Files and tests

| Path | Site / feed | Purpose | Tests |
|------|-------------|---------|-------|
| `ADS-B/healy_repeater/20250623.TSV` | Healy repeater | ADS-B TSV parse and date filtering | `tests/utils/test_adsb.py` |
| `MXAK-AIS-GLBA/MXAK-AIS-GLBA-20250107.csv` | GLBA | Single-day MXAK AIS (e.g. ANNA T) | `tests/utils/ais/test_reader.py`, `test_query.py` |
| `MXAK-AIS-GLBA/MXAK-AIS-GLBA-20240524.csv` | GLBA | Full-day AIS matching NVSPL example date | `tests/utils/ais/test_alignment.py` |
| `MXAK-AIS-GLBA/MXAK-AIS-GLBA-20210512-41.csv` | GLBA | Legacy-formatted MXAK AIS sample | n/a |
| `GLBALSTL_ship_visits.csv` | GLBALSTL | Nearest AIS points at CPA for timing checks | `tests/utils/ais/test_alignment.py` |
| `nvspl_archive/.../NVSPL_GLBALSTL_2024_05_24_*.txt` | GLBALSTL | Hourly SPL for 2024-05-24 | `tests/utils/ais/test_alignment.py` |
| `nvspl_archive/.../NVSPL_DENATRLA_2025_06_23_*.txt` | DENATRLA | Hourly SPL for 2025-06-23 | Offline dev (`DENA_example.config`) |

## Local ground truthing (from repo root)
Running the ground truthing script is a great way to view vessel/flight tracks alongside NVSPL data. The following commands should work out of the box after setting up the project virtual environment. 

**GLBALSTL (AIS):** `-e GLBA_example`

```bash
python -m nps_active_space.scripts.run_ground_truthing \
  -e GLBA_example -u GLBA -s LSTL -y 2024 -t AIS
```

**DENATRLA (ADS-B):** `-e DENA_example`

```bash
python -m nps_active_space.scripts.run_ground_truthing \
  -e DENA_example -u DENA -s TRLA -y 2025 -t ADSB
```

To rebuild `faa/MASTER_DENATRLA_20250623_sample.txt` after changing the example ADS-B day, extract rows from the full FAA `MASTER.txt` for the unique `ICAO_address` values in that day's TSV (see `FAAReleasable` in `utils/models.py`).
