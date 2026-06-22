# Example data

Small, real-format samples for offline unit tests and manual debugging.

**Clock conventions:** MXAK AIS timestamps are UTC-naive; NVSPL and ship-visit times are site-local naive (e.g. GLBALSTL ≈ `America/Juneau`).

| Path | Site / feed | Purpose | Tests |
|------|-------------|---------|-------|
| `ADS-B/healy_repeater/20250623.TSV` | Healy repeater | ADS-B TSV parse and date filtering | `tests/utils/test_adsb.py` |
| `AIS/MXAK-AIS-GLBA-20250107.csv` | GLBA | Single-day MXAK AIS (e.g. ANNA T) | `tests/utils/ais/test_reader.py`, `test_query.py` |
| `AIS/MXAK-AIS-GLBA-20240524.csv` | GLBA | Full-day AIS for clock-alignment regression | `tests/utils/ais/test_alignment.py` |
| `AIS/MXAK-AIS-GLBA-20210512-41.csv` | GLBA | AIS overlapping May 2021 NVSPL window | reserved |
| `GLBALSTL_ship_visits.csv` | GLBALSTL | One-off generated CSV with nearest points from vessels | `tests/utils/ais/test_alignment.py` |
| `NVSPL/GLBALSTL_2024/NVSPL_GLBALSTL_2024_05_24_*.txt` | GLBALSTL | Hourly sound pressure levels for 2024-05-24 | `tests/utils/ais/test_alignment.py` (hours 03, 07, 09) |
