@echo off
setlocal
cd /d "%~dp0.."

if exist ".venv\Scripts\activate.bat" call ".venv\Scripts\activate.bat"
set PYTHONUNBUFFERED=1

REM Mac-parity: example_data site + committed annotations (see AK_example.config).
REM Copy nps_active_space\config\AK_example.config to AK.config and set nmsim/aam paths.

set "OUTDIR=docker"
set "METRICS=python -u docs\tmp_platform_run_metrics.py"
set "OUT=%OUTDIR%\denatrla_1500m_macparity_metrics.json"
set "JOB=--job-file docs\tmp_platform_aam_macparity_job.json --launcher windows-native"
set "GEN=nps_active_space\scripts\generate_active_space.py"
set "BASE=-e AK -u DENA -s TRLA -y 2025 -l 1500 --headings 0 --omni-min 0 --omni-max 0 --density 10 --cleanup --annotation-file DENATRLA2025_saved_annotations.geojson"

echo Mac-parity 1500 m smoke ^(example_data, density 10, omni +000^) — AAM then NMSim
if exist "%OUT%" del "%OUT%"
if exist "%OUTDIR%\denatrla_1500m_macparity_aam.json" del "%OUTDIR%\denatrla_1500m_macparity_aam.json"
if exist "%OUTDIR%\denatrla_1500m_macparity_nmsim.json" del "%OUTDIR%\denatrla_1500m_macparity_nmsim.json"

%METRICS% --out %OUT% --label aam %JOB% -- %GEN% --model aam %BASE% --results-out %OUTDIR%\denatrla_1500m_macparity_aam.json
if errorlevel 1 exit /b 1
%METRICS% --out %OUT% --label nmsim %JOB% -- %GEN% %BASE% --results-out %OUTDIR%\denatrla_1500m_macparity_nmsim.json
exit /b %errorlevel%
