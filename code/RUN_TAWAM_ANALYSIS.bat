@echo off
setlocal
cd /d "%~dp0"

echo ================================================================
echo Tawam Cox / XGBoost-Cox journal reanalysis
echo ================================================================

where python >nul 2>nul
if errorlevel 1 (
  echo ERROR: Python was not found on PATH.
  echo Install Python 3.11 or 3.12 from python.org, select "Add Python to PATH",
  echo then run this file again.
  pause
  exit /b 1
)

if not exist .venv (
  echo Creating virtual environment...
  python -m venv .venv
  if errorlevel 1 goto :fail
)

call .venv\Scripts\activate.bat
if errorlevel 1 goto :fail

python -m pip install --upgrade pip
if errorlevel 1 goto :fail

python -m pip install -r requirements.txt
if errorlevel 1 goto :fail

echo.
echo Running analysis. This can take several minutes because bootstrap CIs use 10,000 resamples.
python run_tawam_journal_reanalysis.py --bootstrap 10000
if errorlevel 1 goto :fail

echo.
echo ================================================================
echo COMPLETE. Please zip the RESULTS folder and upload it to ChatGPT.
echo ================================================================
pause
exit /b 0

:fail
echo.
echo ANALYSIS FAILED. Please copy the full error text or take a screenshot and send it to ChatGPT.
pause
exit /b 1
