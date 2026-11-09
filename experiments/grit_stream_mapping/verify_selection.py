"""Cross-check sample nearest distances without record bbox filtering or STRtree."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import fsspec
import pyarrow.parquet as pq
from analyze import footprint
from pyproj import Transformer
from shapely import from_wkb, transform


def verify(output: Path) -> None:
    receipt = json.loads((output / "river-acquisition.json").read_text())
    cycles = json.loads((output / "inventory.json").read_text())["cycles"]
    project = Transformer.from_crs(4326, 8857, always_xy=True)
    checks = []
    with fsspec.open(
        receipt["url"], "rb", block_size=2**20, headers={"If-Match": receipt["etag"], "Accept-Encoding": "identity"}
    ) as stream:
        parquet = pq.ParquetFile(stream)
        tables = [
            parquet.read_row_group(i, columns=["id", "source_id", "level", "geometry"]) for i in receipt["row_groups"]
        ]
        rows = [r for table in tables for r in table.to_pylist()]
        geoms = [transform(from_wkb(r["geometry"]), project.transform, interleaved=False) for r in rows]
        for i in [3, 36, 100, 293]:
            cell = footprint(*cycles[i]["members"][0]["cell"])
            distances = [cell.distance(g) for g in geoms]
            index = min(range(len(rows)), key=lambda j: distances[j])
            expected = cycles[i]["members"][0]["distance_m"]
            assert distances[index] == expected
            checks.append(
                {
                    "cycle_id": i + 1,
                    "nearest_source_id": rows[index]["source_id"],
                    "nearest_level": rows[index]["level"],
                    "distance_m": distances[index],
                    "inventory_m": expected,
                }
            )
    result = {"unfiltered_rows": len(rows), "row_groups": receipt["row_groups"], "checks": checks}
    (output / "unfiltered-check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    verify(parser.parse_args().output)
