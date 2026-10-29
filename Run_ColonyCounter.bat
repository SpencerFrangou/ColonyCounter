@echo off
rem ColonyCounter: ONE double-click does everything (Windows).
cd /d "%~dp0"
echo == ColonyCounter ==
where python >nul 2>nul || (echo Python not found. Install Python 3.10+ from python.org, tick "Add python.exe to PATH", and run again. & pause & exit /b 1)
if not exist .venv (echo first run: creating the private environment & python -m venv .venv || (echo could not create .venv & pause & exit /b 1))
call .venv\Scripts\activate.bat
python colonycounter\ensure_env.py || (pause & exit /b 1)
for /f "usebackq delims=" %%F in (`powershell -NoProfile -Command "Add-Type -AssemblyName System.Windows.Forms; $d=New-Object System.Windows.Forms.FolderBrowserDialog; $d.Description='Choose the folder that contains the scanned plates'; if($d.ShowDialog() -eq 'OK'){$d.SelectedPath}"`) do set FOLDER=%%F
if "%FOLDER%"=="" (echo no folder chosen & exit /b 0)
echo folder: %FOLDER%
if not exist "%FOLDER%\PlateMap.xlsx" (
  python colonycounter\cfa_platemap.py template "%FOLDER%" && start "" "%FOLDER%\PlateMap.xlsx"
  echo PlateMap.xlsx was created and opened. Fill it in, save, then double-click Run_ColonyCounter.bat again. & pause & exit /b 0
)
python colonycounter\cfa_platemap.py check "%FOLDER%\PlateMap.xlsx" || (echo Fix the plate map and run again. & pause & exit /b 1)
for /f %%T in ('powershell -NoProfile -Command "Get-Date -Format yyyy-MM-dd_HHmm"') do set STAMP=%%T
set OUT=%FOLDER%\Results_%STAMP%
mkdir "%OUT%" 2>nul
python colonycounter\cfa_run.py "%FOLDER%\PlateMap.xlsx" --out "%OUT%" > "%OUT%\run_log.txt" 2>&1
type "%OUT%\run_log.txt" | findstr /v findfont
mkdir "%OUT%\code_used" 2>nul & copy colonycounter\*.py "%OUT%\code_used\" >nul & copy colonycounter\requirements.txt "%OUT%\code_used\" >nul
copy "%FOLDER%\PlateMap.xlsx" "%OUT%\PlateMap_used.xlsx" >nul
start "" "%OUT%\report.html"
echo Done. Results are in %OUT%
pause
