@echo off
setlocal
cd /d "%~dp0.."

if exist ".venv\Scripts\activate.bat" call ".venv\Scripts\activate.bat"
set PYTHONUNBUFFERED=1

REM Quick 3D pipeline smoke on example_data DENATRLA (committed annotations).
REM Copy nps_active_space\config\AK_example.config to AK.config; set [project] nmsim/aam.
REM
REM Usage:
REM   docs\tmp_windows_aam_3d_example.bat          both models (default)
REM   docs\tmp_windows_aam_3d_example.bat aam      AAM only
REM   docs\tmp_windows_aam_3d_example.bat nmsim    NMSim only
REM
REM Macparity mesh/gain per layer (d10, omni +000); 3 layers 1200-1800 m.
REM Full standard job: docs\tmp_windows_aam_3d.bat

set "GEN3D=nps_active_space\scripts\generate_3d_active_space.py"
set "SITE=-e AK -u DENA -s TRLA -y 2025 -a nvspl"
set "ALT=--min-altitude 1200 --max-altitude 1800"
set "ANNOT=DENATRLA2025_saved_annotations.geojson"
set "ANNOT_FULL=example_data\site_projects\DENATRLA\%ANNOT%"
set "EXTRA=--headings 0 --omni-min 0 --omni-max 0 --density 10 --cleanup --annotation-file %ANNOT%"

set "RUN_AAM=1"
set "RUN_NMSIM=1"
if /i "%~1"=="aam" set "RUN_NMSIM=0"
if /i "%~1"=="nmsim" set "RUN_NMSIM=0"
if /i "%~1"=="both" (
  set "RUN_AAM=1"
  set "RUN_NMSIM=1"
)

echo DENATRLA 3D example smoke: layers 1200/1500/1800 m, density 10, omni +000
echo Annotations: %ANNOT_FULL%
echo.

if not exist "%ANNOT_FULL%" (
  echo ERROR: committed example annotations not found: %ANNOT_FULL%
  exit /b 1
)

if "%RUN_AAM%"=="1" (
  echo === AAM 3D example ^(generate + batch + fit^) ===
  python -u %GEN3D% %SITE% --model aam %ALT% %EXTRA%
  if errorlevel 1 exit /b 1
)

if "%RUN_NMSIM%"=="1" (
  echo === NMSim 3D example ^(generate + batch + fit^) ===
  python -u %GEN3D% %SITE% --model nmsim %ALT% %EXTRA%
  if errorlevel 1 exit /b 1
)

echo.
echo Done. Check:
echo   example_data\site_projects\fits.csv
echo   example_data\site_projects\DENATRLA\Output_Data\aam\ACTIVESPACES\DENATRLA2025_*m
echo   example_data\site_projects\DENATRLA\Output_Data\nmsim\ACTIVESPACES\DENATRLA2025_*m
echo Viz:
echo   python nps_active_space\scripts\viz.py DENATRLA2025 -e AK -s -a --compare --terraced --annotation-file %ANNOT%
exit /b 0
