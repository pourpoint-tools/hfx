"""Retrieve only the two official FABDEM members; ZipFile verifies CRC on read."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import requests
from remotezip import RemoteZip

SOURCE = "https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn/"
ARCHIVE = SOURCE + "N30W090-N40W080_FABDEM_V1-2.zip"
MEMBERS = [f"N{lat}W{lon:03d}_FABDEM_V1-2.tif" for lat in (31, 32) for lon in (90, 89)]


def acquire(destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    response = requests.get(SOURCE + "readme.txt", timeout=120)
    response.raise_for_status()
    (destination / "readme.txt").write_bytes(response.content)
    records = []
    with RemoteZip(ARCHIVE) as archive:
        for name in MEMBERS:
            info = archive.getinfo(name)
            data = archive.read(name)
            path = destination / name
            if path.exists() and path.read_bytes() != data:
                raise ValueError(f"Existing input differs: {path}")
            path.write_bytes(data)
            records.append(
                {"member": name, "bytes": len(data), "zip_crc32": info.CRC, "sha256": hashlib.sha256(data).hexdigest()}
            )
    (destination / "provenance.json").write_text(
        json.dumps(
            {
                "version": "FABDEM V1-2",
                "archive": ARCHIVE,
                "readme": SOURCE + "readme.txt",
                "readme_sha256": hashlib.sha256(response.content).hexdigest(),
                "extraction": "HTTP range reads with remotezip; ZIP CRC verified on full member read",
                "members": records,
                "use": "Personal non-commercial testing only; CC BY-NC-SA 4.0; see readme.txt for Copernicus attribution",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    acquire(parser.parse_args().destination)
