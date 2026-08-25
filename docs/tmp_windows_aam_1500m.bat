@echo off
setlocal
cd /d "%~dp0.."

if exist ".venv\Scripts\activate.bat" call ".venv\Scripts\activate.bat"

set "PY=python -u nps_active_space\scripts\generate_active_space.py"
set "BASE=-e windows -u DENA -s TRLA -y 2025 -l 1500 --headings 0 --omni-min 0 --cleanup --annotation-file DENATRLA2025_saved_annotations.geojson"

if /I "%~1"=="smoke" goto smoke

echo Full 1500 m layer ^(density 48, omni 0-2^)
%PY% --model aam %BASE% --omni-max 2 --results-out denatrla_1500m_aam.json
if errorlevel 1 exit /b 1
%PY% %BASE% --omni-max 2 --results-out denatrla_1500m_nmsim.json
exit /b %errorlevel%

:smoke
echo Smoke 1500 m ^(density 10, omni 0^)
%PY% --model aam %BASE% --omni-max 0 --density 10 --results-out denatrla_1500m_aam.json
if errorlevel 1 exit /b 1
%PY% %BASE% --omni-max 0 --density 10 --results-out denatrla_1500m_nmsim.json
exit /b %errorlevel%
