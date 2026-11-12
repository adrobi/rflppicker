#!/usr/bin/env bash
# Build the Linux folder dist/RFLP-Picker with PyInstaller (see requirements-build.txt).
set -euo pipefail
cd "$(dirname "$0")"

python=python3
for candidate in .venv/bin/python .venv313/bin/python; do   # prefer a project environment
  if [ -x "$candidate" ]; then python="$candidate"; break; fi
done

"$python" -m PyInstaller --noconfirm --clean --noupx --onedir --windowed \
  --name RFLP-Picker \
  --collect-all primer3 \
  qtprimer3_vis.py

cp -r names dist/RFLP-Picker/names
# The build contains PyQt6 and Primer3 (GPL): the license text goes with it.
cp LICENSE dist/RFLP-Picker/
# The instructions stay on GitHub; drop the READMEs copied from third-party packages.
find dist/RFLP-Picker -type f -iname 'readme*' -delete

echo "Build complete: dist/RFLP-Picker/RFLP-Picker"
