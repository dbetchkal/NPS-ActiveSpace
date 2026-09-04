@echo off
setlocal
REM ============================================================================
REM AAM site diagnostics — paste output back for debugging run count / timing.
REM
REM Usage:
REM   docs\tmp_aam_site_report.bat
REM   docs\tmp_aam_site_report.bat "T:\path\to\site"
REM ============================================================================

set "SITE=T:\ResMgmt\Sound\NMSim_Partial_Extract\NMSim\01 SITES\DENASUSH"
if not "%~1"=="" set "SITE=%~1"

set "PS1=%~dp0tmp_aam_site_report.ps1"
if not exist "%PS1%" (
  echo ERROR: missing %PS1%
  exit /b 1
)

echo.
powershell -NoProfile -ExecutionPolicy Bypass -File "%PS1%" -Site "%SITE%"
echo.

endlocal
