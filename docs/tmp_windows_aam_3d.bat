@echo off
setlocal
cd /d "%~dp0.."

if exist ".venv\Scripts\activate.bat" call ".venv\Scripts\activate.bat"
set PYTHONUNBUFFERED=1

REM 3D active-space pipeline on example_data DENATRLA 2025 (standard density + omni).
REM Copy nps_active_space\config\AK_example.config to AK.config; set [project] nmsim/aam.
REM
REM Usage:
REM   docs\tmp_windows_aam_3d.bat          both models (default)
REM   docs\tmp_windows_aam_3d.bat aam      AAM only (skip NMSim if layers already exist)
REM   docs\tmp_windows_aam_3d.bat nmsim    NMSim only
REM
REM Pipeline per model: generate_3d_active_space.py (batch + fit_3d_active_space.py).
REM Layers: 900–2100 m (300 m step: 900, 1200, 1500, 1800, 2100).
REM Mesh/gain match docs\tmp_windows_aam_1500m.bat full (not macparity smoke).

set "GEN3D=nps_active_space\scripts\generate_3d_active_space.py"
set "SITE=-e AK -u DENA -s TRLA -y 2025 -a nvspl"
set "ALT=--min-altitude 900 --max-altitude 2100"
set "DENSITY=48"
set "EXTRA=--headings 0 --omni-min 0 --omni-max 2 --density %DENSITY% --cleanup --annotation-file DENATRLA2025_saved_annotations.geojson"

set "RUN_AAM=1"
set "RUN_NMSIM=1"
if /i "%~1"=="aam" set "RUN_NMSIM=0"
if /i "%~1"=="nmsim" set "RUN_NMSIM=0"
if /i "%~1"=="both" (
  set "RUN_AAM=1"
  set "RUN_NMSIM=1"
)

echo DENATRLA 3D: layers 900-2100 m, density %DENSITY%, omni 0-2 ^(5 layers; expect hours at d48 vs macparity ~204s AAM / ~108s NMSim for 1 layer d10 omni0^)
echo Existing NMSim layers: example_data\site_projects\DENATRLA\Output_Data\nmsim\ACTIVESPACES\DENATRLA2025_*m
echo Existing AAM layers:   example_data\site_projects\DENATRLA\Output_Data\aam\ACTIVESPACES\DENATRLA2025_*m
echo.

if "%RUN_AAM%"=="1" (
  echo === AAM 3D ^(generate + batch + fit^) ===
  python -u %GEN3D% %SITE% --model aam %ALT% %EXTRA%
  if errorlevel 1 exit /b 1
)

if "%RUN_NMSIM%"=="1" (
  echo === NMSim 3D ^(generate + batch + fit^) ===
  python -u %GEN3D% %SITE% --model nmsim %ALT% %EXTRA%
  if errorlevel 1 exit /b 1
)

echo.
echo Done. Check:
echo   example_data\site_projects\fits.csv
echo   example_data\site_projects\DENATRLA\Output_Data\aam\PRECISION_RECALL\PrecisionRecallPlot_3d_DENATRLA2025_f1p0.png
echo   example_data\site_projects\DENATRLA\Output_Data\nmsim\PRECISION_RECALL\PrecisionRecallPlot_3d_DENATRLA2025_f1p0.png
echo Viz: python nps_active_space\scripts\viz.py DENATRLA2025 -e AK -s -a --compare --terraced
exit /b 0
