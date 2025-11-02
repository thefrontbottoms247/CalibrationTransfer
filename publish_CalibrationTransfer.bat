@echo off
:: ==========================================================
::  CalibrationTransfer Package – Git + Install Automation
:: ==========================================================
::  This script:
::    1. Initializes git if needed
::    2. Links to GitHub repo
::    3. Commits and pushes cleanly
::    4. Installs package in R via devtools
:: ==========================================================

cd /d "C:\Users\wcc93\Dropbox\(R) RStudio\Calibration Transfer"

echo.
echo ===[ STEP 1: Initialize git repository ]===
git init

echo.
echo ===[ STEP 2: Clean up previous cache (optional) ]===
git rm -r --cached . 2>nul

echo.
echo ===[ STEP 3: Create .gitignore file ]===
(
echo # R package + RStudio
echo .Rhistory
echo .RData
echo .Rproj.user/
echo .Rproj
echo .Rbuildignore

echo # OS junk
echo .DS_Store
echo Thumbs.db
echo desktop.ini

echo # Temporary
echo *.tmp
echo *.log
echo tmp/
echo temp/

echo # Compiled/binary
echo *.o
echo *.so
echo *.dll
echo *.exe

echo # Personal directories
echo C:/Users/wcc93/AppData/
echo C:/Users/wcc93/Documents/
echo C:/Users/wcc93/Downloads/
) > .gitignore

echo.
echo ===[ STEP 4: Stage and commit files ]===
git add .
git commit -m "Initial commit: CalibrationTransfer package"

echo.
echo ===[ STEP 5: Link to GitHub remote ]===
git branch -M main
git remote remove origin 2>nul
git remote add origin https://github.com/thefrontbottoms247/CalibrationTransfer.git

echo.
echo ===[ STEP 6: Push to GitHub (will ask for token) ]===
git push -u origin main

echo.
echo ==========================================================
echo ✅ Push complete — now installing from GitHub in R
echo ==========================================================

echo library(devtools)> "%TEMP%\install_calibtransfer.R"
echo Sys.setenv(GITHUB_PAT = "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")>> "%TEMP%\install_calibtransfer.R"
echo devtools::install_github("thefrontbottoms247/CalibrationTransfer", auth_token = Sys.getenv("GITHUB_PAT"))>> "%TEMP%\install_calibtransfer.R"
echo library(CalibrationTransfer)>> "%TEMP%\install_calibtransfer.R"
echo message("✅ CalibrationTransfer installed successfully!")>> "%TEMP%\install_calibtransfer.R"

Rscript "%TEMP%\install_calibtransfer.R"

echo.
echo ==========================================================
echo ✅ All done! Check GitHub:
echo    https://github.com/thefrontbottoms247/CalibrationTransfer
echo ==========================================================
pause
