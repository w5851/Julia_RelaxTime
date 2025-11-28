@echo off
chcp 65001 >nul
REM 2->2 Scattering Visualization System

REM Resolve repo root (two levels up from scripts/server)
pushd "%~dp0..\.."

echo ========================================
echo    2->2 Scattering Web Visualization
echo ========================================
echo.

REM Check if Julia is installed
where julia >nul 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Julia not found
    echo Download: https://julialang.org/downloads/
    pause
    exit /b 1
)

echo [1/3] Checking Julia dependencies...
julia --project=. -e "using Pkg; Pkg.instantiate()"

echo.
echo [2/3] Starting HTTP server...
echo Note: Press Ctrl+C to stop server
echo.
start "Julia Server" julia scripts/server/server_full.jl

timeout /t 3 /nobreak >nul

echo.
echo [3/3] Opening frontend page...
timeout /t 3 /nobreak >nul
start "" http://localhost:8080

echo.
echo ========================================
echo OK System started!
echo.
echo Server: http://localhost:8080
echo Frontend: web\index.html
echo.
echo Usage:
echo    1. Input parameters in browser
echo    2. Click "Calculate" button
echo    3. View 3D visualization
echo.
echo WARNING: Closing this window will stop the server
echo ========================================
echo.
echo Press any key to exit (server will continue running)...
pause >nul
popd
exit
