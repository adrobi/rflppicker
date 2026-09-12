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
  --icon (Join-Path $PSScriptRoot 'assets\app_icon.ico') `
  --collect-all primer3 `
  --runtime-hook rth_qt6_path.py `
  qtprimer3_vis.py

$distDir = Join-Path $PSScriptRoot 'dist\RFLP-Picker'
$namesSource = Join-Path $PSScriptRoot 'names'
$namesTarget = Join-Path $distDir 'names'
if (-not (Test-Path -LiteralPath $namesSource)) {
  throw "Required names directory is missing: $namesSource"
}
Copy-Item -LiteralPath $namesSource -Destination $namesTarget -Recurse -Force

# The release archive keeps its instructions on GitHub. Remove documentation
# files copied from third-party packages so zipping this folder stays clean.
Get-ChildItem -LiteralPath $distDir -Recurse -File |
  Where-Object { $_.Name -match '(?i)^readme(?:\..*)?$' } |
  Remove-Item -Force

Write-Host "Build complete: $distDir\RFLP-Picker.exe"
