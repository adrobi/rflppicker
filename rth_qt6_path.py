"""Make the bundled Qt runtime win over DLLs installed elsewhere on Windows."""

import os
import sys


if sys.platform == "win32":
    bundle_root = getattr(sys, "_MEIPASS", os.path.dirname(sys.executable))
    qt_bin = os.path.join(bundle_root, "PyQt6", "Qt6", "bin")
    if os.path.isdir(qt_bin):
        os.environ["PATH"] = qt_bin + os.pathsep + os.environ.get("PATH", "")
        try:
            os.add_dll_directory(qt_bin)
        except (AttributeError, OSError):
            pass
