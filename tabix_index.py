"""Read VCF ranges through htslib-compatible TBI or CSI binning indices."""

import gzip
import os
import stat
import struct

from Bio.bgzf import BgzfReader


def index_path(path):
    source_mtime = os.stat(path).st_mtime_ns
    for extension in ('.csi', '.tbi'):
        candidate = str(path) + extension
        try:
            info = os.stat(candidate)
        except FileNotFoundError:
            continue
        if stat.S_ISREG(info.st_mode) and info.st_mtime_ns >= source_mtime:
            return candidate
    return None


def _read(stream, size):
    data = stream.read(size)
    if len(data) != size:
        raise ValueError('Повреждён индекс tabix')
    return data


def _unpack(stream, layout):
    return struct.unpack(layout, _read(stream, struct.calcsize(layout)))


def _read_index(path):
    with gzip.open(path, 'rb') as index:
        magic = _read(index, 4)
        if magic == b'TBI\x01':
            (count,) = _unpack(index, '<I')
            shift, depth = 14, 5
            fmt, _, _, _, _, _, size = _unpack(index, '<7I')
            names = _read(index, size)
        elif magic == b'CSI\x01':
            shift, depth, aux_size = _unpack(index, '<3I')
            auxiliary = _read(index, aux_size)
            if len(auxiliary) < 28:
                raise ValueError('Индекс tabix не содержит имён контигов VCF')
            fmt, _, _, _, _, _, size = struct.unpack_from('<7I', auxiliary)
            names = auxiliary[28:28 + size]
            (count,) = _unpack(index, '<I')
        else:
            raise ValueError('Повреждён индекс tabix')
        if fmt != 2 or not (1 <= shift < 64 and 0 <= depth <= 10 and shift + 3 * depth < 64):
            raise ValueError('Индекс tabix не предназначен для VCF')
        names = names.decode('utf-8').rstrip('\0').split('\0') if names else []
        if len(names) != count:
            raise ValueError('Повреждён индекс tabix')
        references = {}
        for name in names:
            (n_bin,) = _unpack(index, '<I')
            bins = {}
            for _ in range(n_bin):
                if magic == b'TBI\x01':
                    number, n_chunk = _unpack(index, '<2I')
                else:
                    number, _, n_chunk = _unpack(index, '<IQI')
                bins[number] = list(struct.iter_unpack('<2Q', _read(index, 16 * n_chunk)))
            if magic == b'TBI\x01':
                (n_offsets,) = _unpack(index, '<I')
                linear = tuple(value for value, in struct.iter_unpack('<Q', _read(index, 8 * n_offsets)))
            else:
                linear = ()
            references[name] = bins, linear
    return shift, depth, references


def _regions_to_chunks(windows, bins, linear, shift, depth):
    max_end = 1 << (shift + 3 * depth)
    for left, right in windows:
        start, stop = max(left - 1, 0), min(right, max_end)
        if start >= stop:
            continue
        linear_interval = start >> shift
        floor = linear[linear_interval] if linear_interval < len(linear) else 0
        for level in range(depth + 1):
            offset = ((1 << (level * 3)) - 1) // 7
            bits = shift + (depth - level) * 3
            first = offset + (start >> bits)
            last = offset + ((stop - 1) >> bits)
            for bin_number in range(first, last + 1):
                for begin, end in bins.get(bin_number, ()):
                    if end > floor:
                        yield max(begin, floor), end


def vcf_chunks(path, index, select_windows):
    """Yield lines from the selected BGZF chunks; callers still filter by POS."""
    shift, depth, references = _read_index(index)
    ranges = []
    for name, (bins, linear) in references.items():
        ranges.extend(_regions_to_chunks(select_windows(name), bins, linear, shift, depth))
    merged = []
    for start, stop in sorted(ranges):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(stop, merged[-1][1])
        else:
            merged.append([start, stop])
    with BgzfReader(path, 'rb') as stream:
        for start, stop in merged:
            stream.seek(start)
            while stream.tell() < stop:
                line = stream.readline()
                if not line:
                    break
                yield line.decode('utf-8', errors='replace')
