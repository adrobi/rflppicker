$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot

python -m PyInstaller --noconfirm --clean --onedir --windowed `
  --name RFLP-Picker `
  qtprimer3_vis.py

Write-Host "Build complete: $PSScriptRoot\dist\RFLP-Picker\RFLP-Picker.exe"
