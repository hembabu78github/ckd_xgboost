@echo off
setlocal
cd /d "%~dp0"

echo ================================================================
echo Tawam Cox / XGBoost-Cox journal reanalysis - FULL RUN
echo Requires Python 3.12 (64-bit)
echo ================================================================

py -3.12 -c "import sys; print(sys.version)" >nul 2>nul
if errorlevel 1 (
  echo ERROR: Python 3.12 was not found.
  echo Install 64-bit Python 3.12, then run this BAT again.
  echo Check installed versions with: py -0p
  pause
  exit /b 1
)

if exist .venv (
  .venv\Scripts\python.exe -c "import sys; raise SystemExit(0 if sys.version_info[:2]==(3,12) else 1)" >nul 2>nul
  if errorlevel 1 (
    echo Removing old virtual environment created with a different Python version...
    rmdir /s /q .venv
  )
)

if not exist .venv (
  echo Creating Python 3.12 virtual environment...
  py -3.12 -m venv .venv
  if errorlevel 1 goto :fail
)

echo Installing dependencies...
.venv\Scripts\python.exe -m pip install --upgrade pip
if errorlevel 1 goto :fail
.venv\Scripts\python.exe -m pip install --only-binary=:all: -r requirements.txt
if errorlevel 1 goto :fail

echo.
echo Running full analysis. Bootstrap CIs use 10000 resamples.
.venv\Scripts\python.exe run_tawam_journal_reanalysis.py --bootstrap 10000
if errorlevel 1 goto :fail

echo.
echo ================================================================
echo COMPLETE. Zip the RESULTS folder and upload it to ChatGPT.
echo ================================================================
pause
exit /b 0

:fail
echo.
echo ANALYSIS FAILED. Copy the full error text and send it to ChatGPT.
pause
exit /b 1
