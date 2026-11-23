"""Tests for GRIT catchment connectivity completeness metadata preparation and verification."""

from __future__ import annotations

import copy
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import cast

import pytest

from hfx.candidates import Fingerprint


def _fingerprint_bytes(data: bytes) -> Fingerprint:
    size = len(data)
    return {
        "size_bytes": size,
        "sha256": hashlib.sha256(data).hexdigest(),
        "head_sha256": hashlib.sha256(data[:65536]).hexdigest(),
        "tail_sha256": hashlib.sha256(data[max(0, size - 65536) :]).hexdigest(),
    }


def _synthetic_grit_manifest() -> dict[str, object]:
    prefix = "vectors/2026-10-18-final/"
    return {
        "format_version": "0.4.0",
        "fabric_name": "grit",
        "fabric_version": "1.0",
        "crs": "EPSG:4326",
        "has_up_area": False,
        "topology": "dag",
        "bbox": [-180.0, -56.0, 180.0, 84.0],
        "unit_count": 22419857,
        "created_at": "2026-10-13T03:16:01.777967+00:00",
        "adapter_version": "grit-original-vectors-1",
        "region": "global",
        "auxiliary": [
            {
                "schema": "hfx.aux.d8_raster.v2",
                "artifacts": {
                    "flow_dir": {"external": "grit-flow-dir"},
                    "flow_acc": {"external": "grit-flow-acc"},
                },
                "metadata": {"crs": "EPSG:8857", "flow_dir_encoding": "grass", "flow_acc_units": "km2"},
            }
        ],
        "external_artifacts": {
            "grit-flow-dir": {
                "url": "https://basin-delineations-public.upstream.tech/grit/hfx-v0.4.0/aux/d8/flow_dir.tif",
                "size_bytes": 35460747402,
                "sha256": "6d03ccdc092b914be0ae41975fbafebb2b02334ff481a1589897539228ab7c65",
                "head_sha256": "6fdfd1be11525d05e25d90a16558f342497566fac4433f18048f6e771baaf7f0",
                "tail_sha256": "fb8952bac4fa1785df636c7e553cd1dfe5d7327a128cc39060659a9c38605f70",
                "retention": "The GRIT HFX publisher retains this exact public object unchanged for the lifetime of this dependent bundle.",
            },
            "grit-flow-acc": {
                "url": "https://basin-delineations-public.upstream.tech/grit/hfx-v0.4.0/aux/d8/flow_acc.tif",
                "size_bytes": 150995308307,
                "sha256": "01e8b76570a8d98026fdad77107ed5fe719490433d0351c43b62aca32ff0c15d",
                "head_sha256": "8883dd98fa4b3e1aa0650c73564f31bbcbec96c1305075df5c712fed28aafd7a",
                "tail_sha256": "bb4b0c3203ecf773c9cd2cb780c0eb8b272cfeb76f005bf5917ec03a6fcc5d12",
                "retention": "The GRIT HFX publisher retains this exact public object unchanged for the lifetime of this dependent bundle.",
            },
        },
        "catchment_connectivity": [
            {"level": 0, "status": "not_supplied"},
            {"level": 1, "status": "supplied"},
            {"level": 2, "status": "supplied"},
        ],
        "source_domain_metadata": {"rivers": {"topology": "dag"}, "river_nodes": {}},
        "categorical_columns": [
            {
                "artifact": prefix + "catchments.parquet",
                "column": "level",
                "physical_type": "int16",
                "codes": [
                    {"code": 0, "meaning": "component"},
                    {"code": 1, "meaning": "segment local drainage unit"},
                    {"code": 2, "meaning": "reach local drainage unit"},
                ],
                "null_meaning": "not permitted",
            }
        ],
        "representation": "spatial-geoparquet-adjacency",
        "representation_version": 1,
        "families": {
            "catchments": {
                "data": prefix + "catchments.parquet",
                "data_sha256": "a" * 64,
                "lookup": prefix + "catchments-lookup.arrow",
                "directory": prefix + "catchments-lookup-directory.json",
                "count": 22419857,
            },
            "rivers": {
                "data": prefix + "rivers.parquet",
                "data_sha256": "b" * 64,
                "lookup": prefix + "rivers-lookup.arrow",
                "directory": prefix + "rivers-lookup-directory.json",
                "count": 22337300,
            },
            "river_nodes": {
                "data": prefix + "river_nodes.parquet",
                "data_sha256": "c" * 64,
                "lookup": prefix + "river_nodes-lookup.arrow",
                "directory": prefix + "river_nodes-lookup-directory.json",
                "count": 22278858,
            },
        },
        "relationships": {
            "catchments": {
                "data": prefix + "catchments-relationships.arrow",
                "directory": prefix + "catchments-relationships-directory.json",
                "count": 22419857,
                "columns": {"upstream_ids": {"target_domain": "catchments"}},
                "connectivity": [
                    {"level": 0, "status": "not_supplied"},
                    {"level": 1, "status": "supplied"},
                    {"level": 2, "status": "supplied"},
                ],
                "connectivity_source": "manifest.json",
            },
            "rivers": {
                "data": prefix + "rivers-relationships.arrow",
                "directory": prefix + "rivers-relationships-directory.json",
                "count": 22337300,
                "columns": {"upstream_ids": {"target_domain": "rivers"}},
            },
            "river_nodes": {
                "data": prefix + "river_nodes-relationships.arrow",
                "directory": prefix + "river_nodes-relationships-directory.json",
                "count": 22278858,
                "columns": {"outgoing_river_ids": {"target_domain": "rivers"}},
            },
        },
        "source_manifest": prefix + "source-manifest.json",
        "input_manifest_sha256": "b9a38083131dcf1a4c7c4e73ebe9a97dff2a4b7ef91a8df53f8bc38cadfa48db",
        "addressing": {"block": "zero-based Parquet row group", "row": "zero-based row within group"},
        "artifacts": {
            prefix + "catchments.parquet": _fingerprint_bytes(b"catchments"),
            prefix + "rivers.parquet": _fingerprint_bytes(b"rivers"),
            prefix + "source-manifest.json": _fingerprint_bytes(b"source-manifest"),
        },
    }


def test_prepare_and_verify_updates_only_connectivity_and_provenance() -> None:
    from scripts.grit_connectivity_metadata import (
        CONNECTIVITY_PROVENANCE_PATH,
        OCTOBER_MANIFEST_SHA256,
        prepare_metadata_update,
        verify_metadata_update,
    )

    previous = _synthetic_grit_manifest()
    previous_bytes = json.dumps(previous, indent=2).encode()
    prev_sha = hashlib.sha256(previous_bytes).hexdigest()

    prepared = prepare_metadata_update(previous_bytes, expected_previous_sha256=prev_sha)
    assert prepared.preserved_manifest_path == f"provenance/manifest-{prev_sha}.json"
    assert prepared.preserved_manifest_bytes == previous_bytes
    assert prepared.provenance_path == CONNECTIVITY_PROVENANCE_PATH

    updated = cast(dict[str, object], json.loads(prepared.updated_manifest_bytes))
    expected_connectivity = [
        {"level": 0, "status": "not_supplied"},
        {"level": 1, "status": "supplied", "completeness": "complete"},
        {"level": 2, "status": "supplied", "completeness": "complete"},
    ]
    assert updated["catchment_connectivity"] == expected_connectivity
    relationships = cast(dict[str, dict[str, object]], updated["relationships"])
    assert relationships["catchments"]["connectivity"] == expected_connectivity
    assert relationships["rivers"] == cast(dict[str, object], previous["relationships"])["rivers"]
    assert relationships["river_nodes"] == cast(dict[str, object], previous["relationships"])["river_nodes"]

    for key, value in previous.items():
        if key not in {"catchment_connectivity", "relationships", "artifacts"}:
            assert updated[key] == value

    artifacts = cast(dict[str, object], updated["artifacts"])
    prev_artifacts = cast(dict[str, object], previous["artifacts"])
    assert set(artifacts) == set(prev_artifacts) | {prepared.preserved_manifest_path, prepared.provenance_path}
    for key, fp in prev_artifacts.items():
        assert artifacts[key] == fp
    assert artifacts[prepared.preserved_manifest_path] == _fingerprint_bytes(previous_bytes)
    assert artifacts[prepared.provenance_path] == _fingerprint_bytes(prepared.provenance_bytes)

    provenance = cast(dict[str, object], json.loads(prepared.provenance_bytes))
    october = cast(dict[str, object], provenance["october_verification"])
    assert october["public_manifest_sha256"] == OCTOBER_MANIFEST_SHA256
    assert provenance["preserved_previous_manifest"] == {
        "path": prepared.preserved_manifest_path,
        "url": "https://basin-delineations-public.upstream.tech/grit/hfx-v0.4.0/" + prepared.preserved_manifest_path,
        "sha256": prev_sha,
        "size_bytes": len(previous_bytes),
    }

    summary = verify_metadata_update(
        previous_bytes,
        prepared.updated_manifest_bytes,
        prepared.preserved_manifest_bytes,
        prepared.provenance_bytes,
        expected_previous_sha256=prev_sha,
    )
    assert summary["status"] == "passed"
    assert summary["reused_artifact_count"] == len(prev_artifacts)


def test_prepare_rejects_wrong_previous_digest_or_mismatched_connectivity() -> None:
    from scripts.grit_connectivity_metadata import prepare_metadata_update

    previous = _synthetic_grit_manifest()
    previous_bytes = json.dumps(previous).encode()
    prev_sha = hashlib.sha256(previous_bytes).hexdigest()

    with pytest.raises(ValueError, match="digest"):
        prepare_metadata_update(previous_bytes, expected_previous_sha256="0" * 64)

    mutated = copy.deepcopy(previous)
    cast(list[dict[str, object]], mutated["catchment_connectivity"])[1]["completeness"] = "complete"
    mutated_bytes = json.dumps(mutated).encode()
    with pytest.raises(ValueError, match="previous catchment_connectivity"):
        prepare_metadata_update(mutated_bytes, expected_previous_sha256=hashlib.sha256(mutated_bytes).hexdigest())

    mutated = copy.deepcopy(previous)
    cast(dict[str, dict[str, object]], mutated["relationships"])["catchments"]["connectivity"] = []
    mutated_bytes = json.dumps(mutated).encode()
    with pytest.raises(ValueError, match="relationship connectivity"):
        prepare_metadata_update(mutated_bytes, expected_previous_sha256=hashlib.sha256(mutated_bytes).hexdigest())

    mutated = copy.deepcopy(previous)
    cast(dict[str, object], mutated["artifacts"])["provenance/connectivity-completeness-2026-11-23.json"] = (
        _fingerprint_bytes(b"existing")
    )
    mutated_bytes = json.dumps(mutated).encode()
    with pytest.raises(ValueError, match="already exists"):
        prepare_metadata_update(mutated_bytes, expected_previous_sha256=hashlib.sha256(mutated_bytes).hexdigest())
    assert prev_sha != "0" * 64


@pytest.mark.parametrize(
    ("mutation", "pattern"),
    [
        ("mismatched_relationship_connectivity", "relationship connectivity"),
        ("completeness_on_not_supplied", "catchment_connectivity"),
        ("mutated_existing_artifact_hash", "unchanged"),
        ("mutated_external_raster", "top-level"),
        ("mutated_preserved_manifest_bytes", "preserved"),
        ("mismatched_provenance_fingerprint", "fingerprint"),
        ("extra_unexpected_top_level_change", "top-level"),
    ],
)
def test_verify_rejects_mutations(mutation: str, pattern: str) -> None:
    from scripts.grit_connectivity_metadata import prepare_metadata_update, verify_metadata_update

    previous_bytes = json.dumps(_synthetic_grit_manifest(), indent=2).encode()
    prev_sha = hashlib.sha256(previous_bytes).hexdigest()
    prepared = prepare_metadata_update(previous_bytes, expected_previous_sha256=prev_sha)

    updated = json.loads(prepared.updated_manifest_bytes)
    preserved_bytes = prepared.preserved_manifest_bytes
    provenance_bytes = prepared.provenance_bytes

    if mutation == "mismatched_relationship_connectivity":
        updated["relationships"]["catchments"]["connectivity"][1].pop("completeness")
    elif mutation == "completeness_on_not_supplied":
        updated["catchment_connectivity"][0]["completeness"] = "complete"
    elif mutation == "mutated_existing_artifact_hash":
        updated["artifacts"]["vectors/2026-10-18-final/catchments.parquet"]["sha256"] = "0" * 64
    elif mutation == "mutated_external_raster":
        updated["external_artifacts"]["grit-flow-dir"]["sha256"] = "0" * 64
    elif mutation == "mutated_preserved_manifest_bytes":
        preserved_bytes = b"{}"
    elif mutation == "mismatched_provenance_fingerprint":
        provenance_bytes = provenance_bytes + b" "
    elif mutation == "extra_unexpected_top_level_change":
        updated["unit_count"] = 1

    with pytest.raises(ValueError, match=pattern):
        verify_metadata_update(
            previous_bytes,
            json.dumps(updated, indent=2).encode(),
            preserved_bytes,
            provenance_bytes,
            expected_previous_sha256=prev_sha,
        )


def test_cli_prepare_and_verify_roundtrip(tmp_path: Path) -> None:
    previous_bytes = json.dumps(_synthetic_grit_manifest(), indent=2).encode()
    prev_sha = hashlib.sha256(previous_bytes).hexdigest()
    input_manifest = tmp_path / "october-manifest.json"
    input_manifest.write_bytes(previous_bytes)
    candidate_dir = tmp_path / "candidate"

    script = Path(__file__).parents[1] / "scripts" / "grit_connectivity_metadata.py"
    result = subprocess.run(
        [
            sys.executable,
            str(script),
            "prepare",
            "--input-manifest",
            str(input_manifest),
            "--output-dir",
            str(candidate_dir),
            "--expected-previous-sha256",
            prev_sha,
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    payload = json.loads(result.stdout)
    assert payload["status"] == "passed"
    assert (candidate_dir / "manifest.json").is_file()
    assert (candidate_dir / f"provenance/manifest-{prev_sha}.json").read_bytes() == previous_bytes
    assert (candidate_dir / "provenance/connectivity-completeness-2026-11-23.json").is_file()

    verify_result = subprocess.run(
        [
            sys.executable,
            str(script),
            "verify",
            "--previous-manifest",
            str(input_manifest),
            "--candidate-dir",
            str(candidate_dir),
            "--expected-previous-sha256",
            prev_sha,
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    assert json.loads(verify_result.stdout)["status"] == "passed"

    overwrite_attempt = subprocess.run(
        [
            sys.executable,
            str(script),
            "prepare",
            "--input-manifest",
            str(input_manifest),
            "--output-dir",
            str(candidate_dir),
            "--expected-previous-sha256",
            prev_sha,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert overwrite_attempt.returncode != 0


def test_committed_provenance_receipt_consistency() -> None:
    from scripts.grit_connectivity_metadata import (
        OCTOBER_MANIFEST_SHA256,
        OCTOBER_SOURCE_MANIFEST_SHA256,
        UPDATED_CONNECTIVITY,
    )

    receipt_path = Path(__file__).parents[1] / "scripts/provenance/grit-connectivity-completeness.json"
    raw = receipt_path.read_text()
    receipt = cast(dict[str, object], json.loads(raw))
    previous = cast(dict[str, object], receipt["previous_public_manifest"])
    october = cast(dict[str, object], receipt["october_verification_evidence"])
    published = cast(dict[str, object], receipt["published_bundle"])
    live = cast(dict[str, object], receipt["live_public_verification"])

    assert previous["sha256"] == OCTOBER_MANIFEST_SHA256
    assert october["source_manifest_sha256"] == OCTOBER_SOURCE_MANIFEST_SHA256
    assert published["catchment_connectivity"] == UPDATED_CONNECTIVITY
    assert published["manifest_sha256"] == live["manifest_sha256"]
    assert live["status"] == "passed"
    assert live["cache"] is False
    assert live["allow_incomplete"] is False
    rhine = cast(dict[str, object], live["rhine_at_basel_segment"])
    assert rhine["subbasin_count"] == 424
    assert cast(dict[str, object], rhine["coverage"])["status"] == "complete"
    reaches = cast(dict[str, dict[str, object]], live["reach_level_outlets"])
    assert cast(dict[str, object], reaches["netherlands_11_reaches"]["coverage"])["status"] == "complete"
    assert cast(dict[str, object], reaches["australia_2_reaches"]["coverage"])["status"] == "complete"
    comp = cast(dict[str, dict[str, object]], live["component_level"])
    assert comp["allow_incomplete_false"]["error_code"] == "selection.connectivity_unavailable"
    assert comp["allow_incomplete_true"]["error_code"] == "selection.connectivity_unavailable"
    for secret_marker in ("AWS_SECRET", "aws_secret_access_key", "BEGIN PRIVATE KEY"):
        assert secret_marker not in raw
