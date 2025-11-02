@echo off
REM === Build and install the CalibrationTransfer package ===

setlocal
set "R_PATH=C:\Program Files\R\R-4.5.1\bin\R.exe"
set "PACKAGE_DIR=C:\Users\wcc93\Dropbox\(R) RStudio\Calibration Transfer"
set "LOG_FILE=C:\Users\wcc93\Dropbox\(R) RStudio\build_log.txt"

echo.
echo =============================================
echo   Building and installing CalibrationTransfer
echo =============================================

REM Use delayed expansion for proper quotes
"%R_PATH%" CMD BATCH --vanilla ^
  "\"setwd('%PACKAGE_DIR%'); devtools::document(); devtools::install(upgrade='never'); cat('✅ Package rebuilt and installed\\n')\"" ^
  "%LOG_FILE%"

echo.
echo Done. Log written to: %LOG_FILE%
pause
endlocal
