$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot

$pythonExe = $null
$pythonArgs = @()

# Prefer a project virtual environment so the build is reproducible.
foreach ($candidate in @(
    (Join-Path $PSScriptRoot '.venv\Scripts\python.exe'),
    (Join-Path $PSScriptRoot '.venv313\Scripts\python.exe')
  )) {
  if (Test-Path -LiteralPath $candidate) {
    $pythonExe = $candidate
    break
  }
}

if (-not $pythonExe -and (Get-Command py -ErrorAction SilentlyContinue)) {
  & py -3.13 -c "import sys" 2>$null
  if ($LASTEXITCODE -eq 0) {
    $pythonExe = 'py'
    $pythonArgs = @('-3.13')
  }
}

if (-not $pythonExe) {
  $pythonExe = (Get-Command python -ErrorAction Stop).Source
}

& $pythonExe @pythonArgs -m PyInstaller --noconfirm --clean --noupx --onedir --windowed `
  --name RFLP-Picker `
  --runtime-hook rth_qt6_path.py `
  qtprimer3_vis.py

Write-Host "Build complete: $PSScriptRoot\dist\RFLP-Picker\RFLP-Picker.exe"
