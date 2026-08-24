# AAM integration — future work notes

AAM is **not** wired into the ActiveSpace pipeline yet (see [`docker/README.md`](../docker/README.md) smoke test only). The work below is tracked on a **separate branch/PR** from Docker+Wine NMSim support (`feature/nmsim-linux-mac`); that branch only includes AAM **runtime plumbing** (stage script, smoke test) — not pipeline integration.

Status updated Aug 2026 against sibling workspaces:
- `~/dev/aam_translator` — AAM input/output library (terrain, `.inp`, `.POI`, run log).
- `~/dev/nmsim-aam-experiments` — NMSim↔AAM comparison harness, `activespace-experiments/` parity runs, tier-4 reciprocal validation.

---

## What active space generation needs from a propagation model

Active space does **not** model a flyover time history. It repeatedly asks:

> If the aircraft were at source position (x, y, z), what would the **fixed microphone** hear in each 1/3-octave band?

NMSim today: **one fixed receiver** (`.sit`) + **many source points** on one `.trj` → `.tis` with one spectral row per source. Audibility is `any(band > max(ambience, hearing threshold))` for bands 20 Hz–12.5 kHz (`ActiveSpaceGenerator._find_audible_points`). No Lmax, SEL, or time dynamics downstream.

**AAM analogue:** `COMPUTEPOI` with **1 POI** (mic) + **N track points** (each source position), **not** `COMPUTEGRD`. GRD maps ground receivers from one flyover — wrong geometry and wrong output type (scalar event metrics, no per-band spectra at the mic).

See `activespace-experiments/README.md` (NMSim reciprocal two-point ridge) and tier-4 AAM reciprocal runs in `nmsim-aam-experiments` for validation geometry.

---

## `aam_translator` — sibling library (v0.1.0)

Not yet vendored or depended on in this repo. Install: `pip install -e ../aam_translator` (or GitHub URL in that repo's README).

### What it provides today

| Layer | API | Role |
|-------|-----|------|
| **Write** | `write_terrain`, `write_inp`, `write_aam_inputs` | DEM + AOI → `.ELV`/`.IMP` + `COMPUTEPOI` `.inp` |
| **Write** | `hop_speed_kn` | Speed for multi-point tracks (AAM rejects `speed_kn=0` on N>1 tracks) |
| **Read** | `read_poi`, `PoiTimeHistory` | Tecplot `.POI` → time × bands per receiver zone |
| **Read** | `read_run_log`, `AamRunLog` | `READ ERROR`, analysis track table, terrain/impedance extents |
| **Read** | `read_poi_summary_csv` | `.Single.POI.csv` event metrics (Lmax, etc.) |
| **Read** | `read_nmbgf_grid` | `.ELV`/`.IMP` payload round-trip (not AAM `.GRD` output) |
| **Align** | `assert_track_alignment`, `arrival_time_residuals` | Verify `.POI` row `k` ↔ track point `k` |
| **Coords** | `lonlat_to_model_ft`, `read_nmbgf_header` | WGS84 ↔ AAM model feet |

`nmsim-aam-experiments` depends on it for all AAM case generation; it superseded that repo's older `AAM_inp` tooling.

### AAM runtime gotchas (from tier-4 experiments)

| Topic | Detail |
|-------|--------|
| **Track cap** | **400 points max** per `ONE TRACK` (not manual's 500). `MAX_TRACK_POINTS = 400` in `aam_translator`. |
| **NMSim batch size** | ~4000 points/job today → expect **~10× more AAM runs** per mesh batch unless batching changes. |
| **Silent failure** | Exceeding limits or bad decks: AAM exits **0**, writes **no `.POI`**, logs `READ ERROR`. Check `read_run_log(...).ok`, not exit code. |
| **Multi-point speed** | Cannot use `speed_kn=0` with N>1 track points (`INTRTIME`). Use `hop_speed_kn(track, terrain)` — modeling device; omni sources only stay acoustically inert if sphere has single airspeed. |
| **Row order** | `.POI` row `k` = track point `k` by position only. `time_s` is **arrival** time — often non-monotonic for scattered sources. **Never sort by time.** Use `assert_track_alignment` with `DIAGNOSTICS` log. |
| **Band gap** | AAM POI tops ~10 kHz (bands 10–40); audibility uses up to **12.5 kHz** — policy needed in adapter or `_find_audible_points`. |

### Still out of scope for `aam_translator`

- AAM binary execution (Wine/Docker stays here).
- `PropagationModel` adapter and ActiveSpace batching/audibility/mesh logic.
- NetCDF / NPD source mapping (`source_id` is opaque).
- `COMPUTEGRD` writer or `.GRD` reader (different NMBGF dialect).
- Terrain cell-size override (coarsening knob still deferred in library).

Details: `~/dev/aam_translator/docs/activespace_integration_tasks.md`.

---

## Validation status

| Check | Status |
|-------|--------|
| AAM runs under Docker+Wine | ✅ `docker/validate_aam_smoke.py` |
| Omni source refs matched (NMSim ↔ AAM) | ✅ `nmsim-aam-experiments/reports/source_semantics_validation.md` |
| NMSim reciprocal two-point ridge (ActiveSpace writers) | ✅ Tier 2/3 in `activespace-experiments/` (~−37 dB east−west) |
| AAM `COMPUTEPOI` reciprocal two-point (row-aligned) | ✅ Tier 4 in `nmsim-aam-experiments` (`activespace-experiments/runs/tier4*`) |
| Full mesh scale (48×48 lattice, batched AAM) | ❌ Not yet |
| Adapter wired into `ActiveSpaceGenerator` | ✅ `PropagationModel` + `AamPropagationModel` on `feature/aam-propagation-model` |

**Mode decision:** `COMPUTEPOI` is sufficient for core active space; `COMPUTEGRD` is not on the critical path.

---

## Terrain grid resolution (perf) — still deferred

When AAM replaces NMSim, each run loads `.ELV`/`.IMP` terrain grids (separate from `SETUP PARA` calculating Δ).

**Early spike** (single POI, 1000 m AGL, `nmsim-aam-experiments`, Aug 2026):

| Terrain cell | vs baseline wall time | Δ Lmax_dBA |
|--------------|----------------------|------------|
| ~98 ft (30 m DEM) | 1.00× | — |
| ~197 ft (2×) | ~1.0× | −0.03 |
| ~394 ft (4×) | ~1.3× | −0.07 |
| ~787 ft (8×) | ~1.3× | −0.12 |

No re-benchmark at ActiveSpace scale (many sources, lower altitude) yet. `aam_translator` still derives terrain cell size from native DEM posting; coarsening needs a `cell_size_m` override in that library or a prep-step — tracked as deferred there.

Background: `nmsim-aam-experiments/notes/aam_nmbgf.md`, `notes/aam_inp_format.md`.

---

## Other integration threads

- **Source mapping** — Omni validated (~0.8 dB flat residual with matched 1000 ft refs). Real aircraft / NPD-library equivalence open. Implementation reference: `nmsim-aam-experiments/compare/scenario/source.py`.
- **Performance** — Single-POI Docker+Wine runs are ~parity with NMSim (~6–8 s). Real cost is **many batched runs** (400 points/run) for full meshes, not one overflight POI.
- **Terrain CRS** — NMSim `.flt` must be geographic EPSG:4269. Our `setup/elevation.py` already does this (`NMSIM_DST_CRS`). AAM uses `aam_translator` AEQD/NMBGF path — separate from GridFloat.

---

## Next steps (this repo)

Work on branch **`feature/aam-propagation-model`** (or similar), rebased onto `main` after `feature/nmsim-linux-mac` merges.

### Phase 1 — Adapter skeleton — **done on `feature/aam-propagation-model`**

1. ✅ Optional `aam_translator` dependency (`pyproject.toml` `[aam]` / Docker `pip install -e ".[aam]"`).
2. ✅ `PropagationModel` protocol, `NmsimPropagationModel` (extracted), `AamPropagationModel`:
   - `write_terrain` + `write_inp` + `hop_speed_kn`, AAM shim, `read_run_log`, `read_poi`, `assert_track_alignment`.
   - Maps to NMSim TIS-shaped DataFrame (`Xpos`, `Ypos`, `Zpos`, `A`, bands `"10"`…`"12500"`; 12.5 kHz = `NaN`).
3. ✅ Unit tests (`tests/active_space/test_aam_propagation_model.py`); Docker test `docker/validate_aam_propagation_model.py`.

`ActiveSpaceGenerator` accepts optional `propagation_model` (defaults to `NmsimPropagationModel`); batch cap uses `max_points_per_run` (4000 NMSim / 400 AAM).

### Phase 2 — Wire into `ActiveSpaceGenerator` — **partial**

1. ✅ Constructor injection of `AamPropagationModel`; `_run_propagation_model` replaces NMSim-only subprocess path.
2. ✅ Batching at `propagation_model.max_points_per_run` via `_preprocess_source_points`.
3. ✅ Cache CSV, audibility, contour, polygon unchanged.
4. ✅ 12.5 kHz: AAM POI has no band 41 → `NaN` in adapter (audibility uses `max(ambience, hearing threshold)` per band).
5. ❌ CLI/config flag to select AAM in `generate_active_space.py` (still NMSim-only by default).

### Phase 3 — Scale validation and perf

1. Coarse mesh smoke (e.g. `src_pt_density=8`) on DENATRLA with AAM vs NMSim — compare `tested_pts` / active space polygon.
2. Terrain coarsening re-benchmark once adapter exists.
3. Document runtime expectations (runs per active space × headings × gains).

### Explicitly not in initial AAM integration

- `generate_active_space` full gain sweep as default AAM path until Phase 2 is stable.
- `COMPUTEGRD` / ground noise maps / DNL products.
- Replacing NMSim as default propagation model in config until validation passes.

---

## Cross-references

| Resource | Purpose |
|----------|---------|
| `~/dev/aam_translator/README.md` | Install, POI reading, `hop_speed_kn`, alignment |
| `~/dev/aam_translator/docs/activespace_integration_tasks.md` | Library scope, caps, deferred items |
| `~/dev/nmsim-aam-experiments/plans/STATUS_aam_nmsim_comparison.md` | Comparison harness done/open |
| `~/dev/nmsim-aam-experiments/activespace-experiments/README.md` | NMSim ActiveSpace-faithful validation |
| `~/dev/nmsim-aam-experiments/reports/source_semantics_validation.md` | Omni source alignment |
| `~/dev/nmsim-aam-experiments/reports/nmsim_terrain_shielding.md` | GridFloat CRS / shielding |
| [`docker/README.md`](../docker/README.md) | AAM smoke test, `-m aam` |
