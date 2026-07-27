@echo off
REM Launch the MassCube GUI.
REM First run: pip install -r requirements.txt
pushd "%~dp0"
python masscube_gui.py
popd
