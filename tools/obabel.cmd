@echo off
REM Delegate to obabel_cli.py using the active `python` (use the same env as chem-workflow).
python "%~dp0obabel_cli.py" %*
