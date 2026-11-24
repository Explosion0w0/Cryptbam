"""BGZF blocks plus empty "padding" blocks of an exact byte size.

Empty BGZF blocks are legal anywhere, so any subset of padded runs still
concatenates into a valid BAM. htslib requires XLEN == 6, so the padding lives
inside the deflate stream (empty stored / fixed-Huffman blocks), not in EXTRA.
"""
import struct
import zlib

EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
_TAIL = {2: "0300", 3: "020c00", 4: "02083000", 5: "010000ffff", 6: "02040000ffff"}  # inflate to b""
_STORED = bytes.fromhex("000000ffff")                                               # 5-byte empty block


def _block(cdata, payload=b""):
    size = 26 + len(cdata)
    assert size <= 65536
    return (b"\x1f\x8b\x08\x04\0\0\0\0\0\xff\x06\0BC\x02\0" + struct.pack("<H", size - 1) + cdata
            + struct.pack("<II", zlib.crc32(payload), len(payload)))


def compress(data):
    out = []
    for i in range(0, len(data), 0xff00):
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        chunk = data[i:i + 0xff00]
        out.append(_block(c.compress(chunk) + c.flush(), chunk))
    return b"".join(out)


def padding(n):
    """Exactly n bytes (n == 0 or n >= 28) of empty BGZF blocks."""
    out = []
    while n:
        size = min(n, 65536)
        if n - size and n - size < 28:
            size = n - 28
        c = size - 26                      # deflate bytes needed
        r = (c - 2) % 5 + 2
        out.append(_block(_STORED * ((c - r) // 5) + bytes.fromhex(_TAIL[r])))
        n -= size
    return b"".join(out)


def pad_run(content, seg, filler=padding, min_gap=28):
    """1 junk byte + content + filler, total a multiple of seg."""
    n = -(-(len(content) + 1) // seg) * seg
    while 0 < n - len(content) - 1 < min_gap:
        n += seg
    return b"\0" + content + filler(n - len(content) - 1)
