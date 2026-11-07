@echo off
rem Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
rem Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
rem for valuable scientific input and feedback during the development of CHASE.
rem Start CHASE from its conda environment. Run from an Anaconda/Miniforge Prompt, or double-click.
cd /d "%~dp0"
where conda >nul 2>nul
if errorlevel 1 (
  for %%D in ("%USERPROFILE%\miniforge3" "%USERPROFILE%\mambaforge" "%USERPROFILE%\miniconda3" "%USERPROFILE%\anaconda3") do (
    if exist "%%~D\condabin\conda.bat" set "PATH=%%~D\condabin;%PATH%"
  )
)
where conda >nul 2>nul
if errorlevel 1 (
  echo conda was not found. Install Miniforge from https://conda-forge.org/download/ and try again.
  pause
  exit /b 1
)
call conda env list | findstr /b /c:"chase " >nul
if errorlevel 1 (
  echo Creating the 'chase' environment ^(first run only^)...
  call conda env create -f environment.yml
)
call conda activate chase
python -m chase %*
pause
