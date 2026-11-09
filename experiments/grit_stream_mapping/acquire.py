"""Select public river records by native-cell envelope, retaining complete lines."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import fsspec
import pyarrow as pa
import pyarrow.parquet as pq
import requests
from pyproj import Transformer

BASE = "https://basin-delineations-public.upstream.tech/grit/hfx-v0.4.0/"


def acquire(output: Path, bounds: tuple[float, float, float, float]) -> None:
    output.mkdir(parents=True, exist_ok=True)
    response = requests.get(BASE + "manifest.json")
    response.raise_for_status()
    (output / "manifest.json").write_bytes(response.content)
    manifest = response.json()
    path = manifest["families"]["rivers"]["data"]
    identity = manifest["artifacts"][path]
    url = BASE + path
    headers = {"Accept-Encoding": "identity"}
    head = requests.get(url, headers={**headers, "Range": "bytes=0-65535"})
    head.raise_for_status()
    etag = head.headers["ETag"]
    tail = requests.get(url, headers={**headers, "Range": f"bytes={identity['size_bytes'] - 65536}-", "If-Match": etag})
    assert head.status_code == tail.status_code == 206
    assert hashlib.sha256(head.content).hexdigest() == identity["head_sha256"]
    assert hashlib.sha256(tail.content).hexdigest() == identity["tail_sha256"]
    west, south, east, north = Transformer.from_crs(8857, 4326, always_xy=True).transform_bounds(
        *bounds, densify_pts=101
    )
    selected = []
    groups = []
    with fsspec.open(url, "rb", block_size=2**20, headers={**headers, "If-Match": etag}) as stream:
        parquet = pq.ParquetFile(stream)
        names = [parquet.metadata.row_group(0).column(i).path_in_schema for i in range(parquet.metadata.num_columns)]
        for i in range(parquet.num_row_groups):
            group = parquet.metadata.row_group(i)
            stats = {
                name: group.column(names.index("bbox." + name)).statistics for name in ("xmin", "ymin", "xmax", "ymax")
            }
            if any(stat is None for stat in stats.values()):
                keep = True
            else:
                keep = not (
                    stats["xmin"].min > east
                    or stats["xmax"].max < west
                    or stats["ymin"].min > north
                    or stats["ymax"].max < south
                )
            if keep:
                table = parquet.read_row_group(
                    i,
                    columns=["id", "level", "source_id", "geometry", "bbox", "upstream_node_id", "downstream_node_id"],
                )
                rows = [
                    row
                    for row in table.to_pylist()
                    if row["bbox"]["xmin"] <= east
                    and row["bbox"]["xmax"] >= west
                    and row["bbox"]["ymin"] <= north
                    and row["bbox"]["ymax"] >= south
                ]
                if rows:
                    selected.append(pa.Table.from_pylist(rows, schema=table.schema))
                groups.append(i)
                print(i, len(rows), flush=True)
    pq.write_table(pa.concat_tables(selected), output / "rivers.parquet")
    receipt = {
        "url": url,
        "identity": identity,
        "etag": etag,
        "bounded_fingerprints_verified": True,
        "full_remote_sha256_verified": False,
        "query_bounds_epsg8857": bounds,
        "query_bounds_epsg4326": [west, south, east, north],
        "row_groups": groups,
        "rows": sum(t.num_rows for t in selected),
    }
    (output / "river-acquisition.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--bounds", type=float, nargs=4, required=True)
    args = parser.parse_args()
    acquire(args.output, tuple(args.bounds))
