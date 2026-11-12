"""The description of a GitHub release: the CHANGELOG section of the version and the archives.

Used by .github/workflows/release.yml::

    python tools/release_notes.py v1.5.0 RFLP-Picker-v1.5.0-windows-x64.zip ... > notes.md

The tag must equal the version of the program, otherwise the release is not made.
"""
import hashlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPOSITORY = 'https://github.com/adrobi/rflppicker'


def version_of_tag(tag: str) -> str:
    """'1.5.0' for the tag 'v1.5.0'; exits when the program has another version."""
    source = (ROOT / 'rflp_core.py').read_text(encoding='utf-8')
    version = re.search(r'^APP_VERSION = "([^"]+)"', source, re.MULTILINE).group(1)
    if tag != 'v' + version:
        sys.exit(f'The tag {tag} does not match APP_VERSION {version} in rflp_core.py')
    return version


def section(changelog: str, version: str) -> str:
    """The text under the heading "## <version>" up to the next version heading."""
    match = re.search(rf'^## {re.escape(version)}[ \t]*\n(.*?)(?=^## |\Z)', changelog, re.MULTILINE | re.DOTALL)
    if match is None or not re.fullmatch(r'\d+\.\d+\.\d+', version):
        sys.exit(f'CHANGELOG.md has no section "## {version}"')
    return match.group(1).strip()


def notes(changelog: str, version: str, checksums: dict) -> str:
    lines = [section(changelog, version), '', '## Загрузка · Download', '']
    lines += [f'- `{name}` — SHA-256 `{digest}`' for name, digest in checksums.items()]
    lines += ['',
              'Распакуйте архив целиком и запустите `RFLP-Picker`; папки `_internal` и `names` должны оставаться '
              'рядом. Python не требуется. Референсный FASTA в архив не входит.',
              '',
              'Unpack the whole archive and run `RFLP-Picker`; the `_internal` and `names` folders must stay next '
              'to it. Python is not required. The reference genome FASTA is not included.',
              '',
              f'[README]({REPOSITORY}/blob/v{version}/README.md) · '
              f'[README на русском]({REPOSITORY}/blob/v{version}/README.ru.md) · '
              f'[Изменения]({REPOSITORY}/blob/v{version}/CHANGELOG.md)']
    return '\n'.join(lines) + '\n'


def main(argv) -> int:
    tag, *archives = argv
    version = version_of_tag(tag)
    checksums = {Path(name).name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in archives}
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stdout.write(notes((ROOT / 'CHANGELOG.md').read_text(encoding='utf-8'), version, checksums))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
