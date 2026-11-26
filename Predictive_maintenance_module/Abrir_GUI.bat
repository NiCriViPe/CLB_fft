@echo off
REM Doble clic para abrir el Modulo de Mantenimiento Predictivo
cd /d "%~dp0"
python gui.py
if errorlevel 1 pause
