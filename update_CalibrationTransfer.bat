@echo off
:: ==========================================================
::  CalibrationTransfer One-Click Update Script
:: ==========================================================
::  This script:
::    1. Builds docs (roxygen2)
::    2. Reinstalls locally
::    3. Commits any changes
::    4. Pushes to GitHub
:: ==========================================================

cd /d "C:\Users\wcc93\Dropbox\(R) RStudio\Calibration Transfer"

echo.
echo ===[ STEP 1: Updating package documentation and build ]===
Rscript -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools', repos='https://cloud.r-project.org'); \
  devtools::document(); devtools::install(upgrade='never'); cat('\n✅ Build and install complete.\n')"

echo.
echo ===[ STEP 2: Stage changes for commit ]===
git add .

echo.
echo ===[ STEP 3: Commit changes with timestamp ]===
for /f "tokens=1-3 delims=/ " %%a in ("%date%") do set datestr=%%a-%%b-%%c
git commit -m "Automated update %datestr%"

echo.
echo ===[ STEP 4: Push to GitHub ]===
git push origin main

echo.
echo ==========================================================
echo ✅ CalibrationTransfer package updated and pushed.
echo    https://github.com/thefrontbottoms247/CalibrationTransfer
echo ==========================================================
pause
