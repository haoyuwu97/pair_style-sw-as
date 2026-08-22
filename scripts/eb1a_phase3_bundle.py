from __future__ import annotations

import base64
import hashlib
import re
import runpy
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "audit_inputs" / "phase3_bundle"

SCRIPT_PARTS = [
    "phase3_script.part00.b64",
    "phase3_script.part01.b64",
]
SEED_PARTS = [
    "phase3_seed_v2.part00.b64",
    "phase3_seed_v2.part01.b64",
    "phase3_seed_v2.part02.b64",
    "phase3_seed_v2.part03.b64",
    "phase3_seed_v2.part04.b64",
    "phase3_seed_v2.part05.b64",
    "phase3_seed_v2.part06.b64",
    "phase3_seed_v2.part07.b64",
]

EXPECTED_SCRIPT_SHA256 = "ee0abf4ea08dcda082f6fa146517c6d0dacf334701e20c11d725fce2c9b4948f"
EXPECTED_SEED_SHA256 = "2017fdd9404ab688b5f1794d842cdf1d9c988d39abfeb7ba9dbe8bc667a6666c"
EXPECTED_SEED_B64_SHA256 = "08971407358013542566c35d3ed841ec01f3215c7b177d0eeb5e0ed2563fbd84"
EXPECTED_SCRIPT_BYTES = 39058
EXPECTED_SEED_BYTES = 91326


def read_encoded(parts: list[str]) -> str:
    encoded = "".join((BUNDLE / part).read_text(encoding="utf-8") for part in parts)
    return re.sub(r"\s+", "", encoded)


def decode_bundle(parts: list[str], expected_b64_sha256: str | None = None) -> bytes:
    encoded = read_encoded(parts)
    encoded_sha256 = hashlib.sha256(encoded.encode("ascii")).hexdigest()
    print(
        f"Read {len(parts)} bundle parts: {len(encoded)} base64 characters, "
        f"sha256={encoded_sha256}"
    )
    if expected_b64_sha256 and encoded_sha256 != expected_b64_sha256:
        raise RuntimeError(
            "Encoded bundle SHA-256 mismatch: "
            f"expected {expected_b64_sha256}, got {encoded_sha256}"
        )
    return zlib.decompress(base64.b64decode(encoded, validate=True))


def verify(data: bytes, expected_bytes: int, expected_sha256: str, label: str) -> None:
    actual_sha256 = hashlib.sha256(data).hexdigest()
    if len(data) != expected_bytes:
        raise RuntimeError(
            f"{label} byte-length mismatch: expected {expected_bytes}, got {len(data)}"
        )
    if actual_sha256 != expected_sha256:
        raise RuntimeError(
            f"{label} SHA-256 mismatch: expected {expected_sha256}, got {actual_sha256}"
        )
    print(f"Verified {label}: {len(data)} bytes, sha256={actual_sha256}")


def main() -> None:
    script_data = decode_bundle(SCRIPT_PARTS)
    seed_data = decode_bundle(SEED_PARTS, EXPECTED_SEED_B64_SHA256)

    verify(script_data, EXPECTED_SCRIPT_BYTES, EXPECTED_SCRIPT_SHA256, "Phase 3 script")
    verify(seed_data, EXPECTED_SEED_BYTES, EXPECTED_SEED_SHA256, "Phase 3 seed")

    input_dir = ROOT / "audit_inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    seed_path = input_dir / "phase3_seed.csv"
    generated_script = ROOT / "scripts" / "_eb1a_phase3_generated.py"
    seed_path.write_bytes(seed_data)
    generated_script.write_bytes(script_data)

    print(f"Materialized seed at {seed_path}")
    print(f"Materialized audit script at {generated_script}")
    runpy.run_path(str(generated_script), run_name="__main__")


if __name__ == "__main__":
    main()
