@echo off
setlocal
cd /d "%~dp0.."

if exist ".venv\Scripts\activate.bat" call ".venv\Scripts\activate.bat"
set PYTHONUNBUFFERED=1

REM DENASUSH 2026 smoke (same flags as macparity DENATRLA, different site).
REM Requires example_data\site_projects\DENASUSH\ (or [project] dir) with
REM DENASUSH2026_saved_annotations.geojson, elevation, and site metadata.
REM Copy AK_example.config to AK.config and set nmsim/aam paths + dem if needed.

set "OUTDIR=docker"
set "METRICS=python -u docs\tmp_platform_run_metrics.py"
set "OUT=%OUTDIR%\denasush2026_test_metrics.json"
set "JOB=--job-file docs\tmp_platform_test_denasush2026_job.json --launcher windows-native"
set "GEN=nps_active_space\scripts\generate_active_space.py"
set "BASE=-e AK -u DENA -s SUSH -y 2026 -l 1500 --headings 0 --omni-min 0 --omni-max 0 --density 10 --cleanup --annotation-file DENASUSH2026_saved_annotations.geojson"

echo DENASUSH2026 1500 m smoke ^(density 10, omni +000^) — AAM then NMSim
if exist "%OUT%" del "%OUT%"
if exist "%OUTDIR%\denasush2026_test_aam.json" del "%OUTDIR%\denasush2026_test_aam.json"
if exist "%OUTDIR%\denasush2026_test_nmsim.json" del "%OUTDIR%\denasush2026_test_nmsim.json"

%METRICS% --out %OUT% --label aam %JOB% -- %GEN% --model aam %BASE% --results-out %OUTDIR%\denasush2026_test_aam.json
if errorlevel 1 exit /b 1
%METRICS% --out %OUT% --label nmsim %JOB% -- %GEN% %BASE% --results-out %OUTDIR%\denasush2026_test_nmsim.json
exit /b %errorlevel%
