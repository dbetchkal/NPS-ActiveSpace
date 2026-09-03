@echo off
setlocal
cd /d "%~dp0.."

if exist ".venv\Scripts\activate.bat" call ".venv\Scripts\activate.bat"
set PYTHONUNBUFFERED=1

REM ============================================================================
REM DENASUSH2026 3D pipeline — edit CONFIG below.
REM Pipeline per model: generate_3d_active_space.py (batch + fit_3d_active_space.py).
REM Layers: MIN_ALTITUDE–MAX_ALTITUDE m in 300 m steps (see generate_3d_active_space.py).
REM
REM Usage:
REM   docs\tmp_test_denasush2026.bat          both models (default)
REM   docs\tmp_test_denasush2026.bat aam      AAM only
REM   docs\tmp_test_denasush2026.bat nmsim    NMSim only
REM ============================================================================

REM --- CONFIG (edit these) ---
set "ENV=AK"
set "UNIT=DENA"
set "SITE=SUSH"
set "YEAR=2026"
set "MIN_ALTITUDE=1200"
set "MAX_ALTITUDE=2400"
set "DENSITY=10"
REM Omni ladder step is 0.5 dB (0-2 = gains 0, 0.5, 1, 1.5, 2). Layers with gain 0 only still run for missing gains.
set "OMNI_MIN=0"
set "OMNI_MAX=2"
set "HEADINGS=0"
set "ANNOT=DENASUSH2026_saved_annotations.geojson"
REM Optional: set USE_METRICS=1 to append wall/cpu to docker\denasush2026_test_metrics.json
set "USE_METRICS=0"
REM --- end CONFIG ---

set "GEN3D=nps_active_space\scripts\generate_3d_active_space.py"
set "SITE=-e %ENV% -u %UNIT% -s %SITE% -y %YEAR% -a nvspl"
set "ALT=--min-altitude %MIN_ALTITUDE% --max-altitude %MAX_ALTITUDE%"
set "EXTRA=--headings %HEADINGS% --omni-min %OMNI_MIN% --omni-max %OMNI_MAX% --density %DENSITY% --cleanup --annotation-file %ANNOT%"

set "OUTDIR=docker"
set "METRICS=python -u docs\tmp_platform_run_metrics.py"
set "OUT=%OUTDIR%\denasush2026_test_metrics.json"
set "JOB=--job-file docs\tmp_platform_test_denasush2026_job.json --launcher windows-native"

set "RUN_AAM=1"
set "RUN_NMSIM=1"
if /i "%~1"=="aam" set "RUN_NMSIM=0"
if /i "%~1"=="nmsim" set "RUN_AAM=0"
if /i "%~1"=="both" (
  set "RUN_AAM=1"
  set "RUN_NMSIM=1"
)

echo DENASUSH2026 3D: layers %MIN_ALTITUDE%-%MAX_ALTITUDE% m ^(300 m step^), density %DENSITY%, omni %OMNI_MIN%-%OMNI_MAX%
echo Models: AAM=%RUN_AAM% NMSim=%RUN_NMSIM%
echo Site dir: see [project] dir in %ENV%.config ^(DENASUSH^)
echo.

if "%USE_METRICS%"=="1" if exist "%OUT%" del "%OUT%"

if "%RUN_AAM%"=="1" (
  echo === AAM 3D ^(generate + batch + fit^) ===
  if "%USE_METRICS%"=="1" (
    %METRICS% --out %OUT% --label aam_3d %JOB% -- python -u %GEN3D% %SITE% --model aam %ALT% %EXTRA%
  ) else (
    python -u %GEN3D% %SITE% --model aam %ALT% %EXTRA%
  )
  if errorlevel 1 exit /b 1
)

if "%RUN_NMSIM%"=="1" (
  echo === NMSim 3D ^(generate + batch + fit^) ===
  if "%USE_METRICS%"=="1" (
    %METRICS% --out %OUT% --label nmsim_3d %JOB% -- python -u %GEN3D% %SITE% --model nmsim %ALT% %EXTRA%
  ) else (
    python -u %GEN3D% %SITE% --model nmsim %ALT% %EXTRA%
  )
  if errorlevel 1 exit /b 1
)

echo.
echo Done. Check:
echo   example_data\site_projects\fits.csv  ^(or your [project] dir^)
echo   DENASUSH\Output_Data\aam\ACTIVESPACES\DENASUSH2026_*m
echo   DENASUSH\Output_Data\nmsim\ACTIVESPACES\DENASUSH2026_*m
echo   DENASUSH\Output_Data\{aam,nmsim}\PRECISION_RECALL\PrecisionRecallPlot_3d_DENASUSH2026_f1p0.png
if "%USE_METRICS%"=="1" echo   Metrics: %OUT%
exit /b 0
