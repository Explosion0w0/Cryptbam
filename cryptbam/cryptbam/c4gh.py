"""Minimal Crypt4GH v1 (GA4GH) with several data keys and per-reader edit lists.

File   = "crypt4gh" | version u32=1 | n_packets u32 | packets | segments
Packet = len u32 | method u32=0 | writer_pk 32 | nonce 12 | ChaCha20-Poly1305(payload)
         payload 0: data key  (u32 type=0, u32 method=0, key 32)
         payload 1: edit list (u32 type=1, u32 n, n x u64 skip/keep lengths)
Segment= nonce 12 | ChaCha20-Poly1305(<= 64 KiB plaintext)

Each reader can only open the header packets encrypted to its public key, so it
only learns its own data keys and its own edit list. Byte-compatible with the
reference `crypt4gh` CLI.
"""
import base64
import hashlib
import os
import struct

from cryptography.exceptions import InvalidTag
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey, X25519PublicKey
from cryptography.hazmat.primitives.ciphers.aead import ChaCha20Poly1305
from cryptography.hazmat.primitives.serialization import Encoding, NoEncryption, PrivateFormat, PublicFormat

SEG = 65536                 # plaintext segment size
CSEG = SEG + 28             # nonce + MAC
KEY_MAGIC = b"c4gh-v1"


# ---- key files (Crypt4GH format, unencrypted private keys) ----------------
def _pem(path, kind, body=None):
    if body is None:
        lines = [l.strip() for l in open(path) if l.strip()]
        return base64.b64decode("".join(lines[1:-1]))
    with open(path, "w") as f:
        f.write(f"-----BEGIN CRYPT4GH {kind} KEY-----\n{base64.b64encode(body).decode()}\n"
                f"-----END CRYPT4GH {kind} KEY-----\n")


def _s(b):
    return struct.pack(">H", len(b)) + b


def keygen(sk_path, pk_path):
    sk = X25519PrivateKey.generate()
    raw = sk.private_bytes(Encoding.Raw, PrivateFormat.Raw, NoEncryption())
    _pem(sk_path, "PRIVATE", KEY_MAGIC + _s(b"none") + _s(b"none") + _s(raw))
    os.chmod(sk_path, 0o600)
    _pem(pk_path, "PUBLIC", sk.public_key().public_bytes(Encoding.Raw, PublicFormat.Raw))


def load_pk(path):
    return _pem(path, "PUBLIC")


def load_sk(path):
    d = _pem(path, "PRIVATE")
    if not d.startswith(KEY_MAGIC + _s(b"none") + _s(b"none")):
        raise ValueError(f"{path}: only unencrypted keys are supported (crypt4gh-keygen --nocrypt)")
    return d[len(KEY_MAGIC) + 14:][:32]


# ---- header packets (X25519 + libsodium crypto_kx + ChaCha20-Poly1305) -----
def _kx(dh, reader_pk, writer_pk):
    return hashlib.blake2b(dh + reader_pk + writer_pk, digest_size=64).digest()[:32]


def _pub(sk):
    return sk.public_key().public_bytes(Encoding.Raw, PublicFormat.Raw)


def seal_packet(payload, writer_sk, reader_pk):
    wpk, nonce = _pub(writer_sk), os.urandom(12)
    key = _kx(writer_sk.exchange(X25519PublicKey.from_public_bytes(reader_pk)), reader_pk, wpk)
    body = struct.pack("<I", 0) + wpk + nonce + ChaCha20Poly1305(key).encrypt(nonce, payload, None)
    return struct.pack("<I", len(body) + 4) + body


def open_packet(body, sk_raw):
    if struct.unpack("<I", body[:4])[0] != 0:
        return None
    sk = X25519PrivateKey.from_private_bytes(sk_raw)
    wpk = body[4:36]
    key = _kx(sk.exchange(X25519PublicKey.from_public_bytes(wpk)), _pub(sk), wpk)
    try:
        return ChaCha20Poly1305(key).decrypt(body[36:48], body[48:], None)
    except InvalidTag:
        return None          # packet addressed to someone else


# ---- writer -----------------------------------------------------------------
def edit_list(runs, readable):
    """Skip/keep lengths for one reader. runs = [(key_id, bytes)].
    Every run starts with 1 junk byte that everybody skips: the reference
    crypt4gh (<=1.8.6) wrongly decrypts the last segment of a skip whose length
    is an exact multiple of 64 KiB, so a skip must end inside a readable segment."""
    out, skip = [], 0
    for kid, data in runs:
        if kid in readable:
            out += [skip + 1, len(data) - 1]
            skip = 0
        else:
            skip += len(data)
    assert skip == 0, "last run must be readable by everyone"
    return out


def write(f, runs, keys, recipients):
    """runs: [(key_id, bytes)], each a multiple of SEG except the last.
    keys: key_id -> 32 bytes. recipients: [(public_key, set_of_key_ids)]."""
    wsk = X25519PrivateKey.generate()          # ephemeral writer key
    packets = []
    for pk, kids in recipients:
        packets += [seal_packet(struct.pack("<II", 0, 0) + keys[k], wsk, pk) for k in sorted(kids)]
        el = edit_list(runs, kids)
        payload = struct.pack(f"<II{len(el)}Q", 1, len(el), *el)
        packets.append(seal_packet(payload, wsk, pk))
    f.write(b"crypt4gh" + struct.pack("<II", 1, len(packets)) + b"".join(packets))
    for kid, data in runs:
        aead = ChaCha20Poly1305(keys[kid])
        for i in range(0, len(data), SEG):
            nonce = os.urandom(12)
            f.write(nonce + aead.encrypt(nonce, data[i:i + SEG], None))
    return {"header_packets": len(packets), "segments": sum(-(-len(d) // SEG) for _, d in runs)}


# ---- reader -----------------------------------------------------------------
def read_header(f, sk_raw):
    magic, ver, n = struct.unpack("<8sII", f.read(16))
    if magic != b"crypt4gh" or ver != 1:
        raise ValueError("not a Crypt4GH v1 file")
    keys, edits = [], None
    for _ in range(n):
        ln = struct.unpack("<I", f.read(4))[0]
        p = open_packet(f.read(ln - 4), sk_raw)
        if p is None:
            continue
        typ, x = struct.unpack("<II", p[:8])
        if typ == 0:
            keys.append(ChaCha20Poly1305(p[8:40]))
        elif edits is None:
            edits = list(struct.unpack(f"<{x}Q", p[8:8 + 8 * x]))
        else:
            raise ValueError("more than one edit list for this reader")
    if not keys:
        raise PermissionError("no data key in this file is addressed to this private key")
    return keys, edits


def decrypt(f, out, sk_raw, use_edit_list=True):
    """Decrypt what this key may read. Segments fully skipped by the edit list
    are never touched. With use_edit_list=False every segment is tried (audit)."""
    keys, edits = read_header(f, sk_raw)
    keep, pos = [], 0                          # absolute plaintext keep ranges
    for i, n in enumerate(edits if use_edit_list and edits else []):
        if i % 2:
            keep.append((pos, pos + n))
        pos += n
    if not use_edit_list or not edits:
        keep = [(0, float("inf"))]
    elif len(edits) % 2:
        keep.append((pos, float("inf")))
    st = {"decrypted": 0, "skipped": 0, "failed": 0}
    i = 0
    while cseg := f.read(CSEG):
        a, b = i * SEG, i * SEG + len(cseg) - 28
        i += 1
        want = [(max(s, a) - a, min(e, b) - a) for s, e in keep if s < b and e > a]
        if not want:
            st["skipped"] += 1
            continue
        plain = None
        for k in keys:
            try:
                plain = k.decrypt(cseg[:12], cseg[12:], None)
                break
            except InvalidTag:
                pass
        if plain is None:
            if use_edit_list:
                raise ValueError(f"segment {i - 1}: authentication failed (wrong key or tampered)")
            st["failed"] += 1
            continue
        st["decrypted"] += 1
        for s, e in want:
            out.write(plain[s:e])
    return st
