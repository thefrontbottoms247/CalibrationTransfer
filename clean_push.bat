:: ==========================================================
:: CalibrationTransfer Git Cleanup + Push Script
:: ==========================================================
:: Run this from: C:\Users\wcc93\Dropbox\(R) RStudio\Calibration Transfer
:: Purpose: Clean repo, ignore junk, and push clean package to GitHub
:: ==========================================================

:: Ensure Git works
git --version

:: Stop tracking everything temporarily
git rm -r --cached .

:: Create a proper .gitignore file
(
echo # R package and RStudio ignores
echo .Rhistory
echo .RData
echo .Rproj.user/
echo .Rproj
echo .Rbuildignore

echo.
echo # OS junk
echo .DS_Store
echo Thumbs.db
echo desktop.ini

echo.
echo # Compiled / binary
echo *.o
echo *.so
echo *.dll
echo *.exe
echo *.cpp~
echo *.obj
echo *.lib
echo *.log

echo.
echo # Local and temp
echo data-raw/
echo inst/doc/
echo vignettes/
echo tmp/
echo temp/
echo build/

echo.
echo # Personal system folders
echo C:/Users/wcc93/AppData/
echo C:/Users/wcc93/Documents/
echo C:/Users/wcc93/Downloads/
echo C:/Users/wcc93/OneDrive/
) > .gitignore

:: Stage only clean project files
git add .gitignore
git add DESCRIPTION NAMESPACE README.md
git add R

:: Commit cleanup
git commit -m "Clean repo: added .gitignore and removed AppData artifacts"

:: Fix line-ending warnings permanently
git config core.autocrlf true

:: Push clean version (force overwrite remote)
git push -f origin main

:: Done
echo ===========================================================
echo ✅ Repository cleaned and pushed successfully.
echo Check https://github.com/wcc93/CalibrationTransfer
echo ===========================================================
pause
