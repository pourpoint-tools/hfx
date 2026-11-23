"""Prepare and verify connectivity completeness metadata for the public GRIT HFX 0.4.0 bundle."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import cast

from hfx.candidates import Fingerprint

PUBLIC_BUNDLE_URL = "https://basin-delineations-public.upstream.tech/grit/hfx-v0.4.0/"
OCTOBER_MANIFEST_SHA256 = "04bbfb50907feaf339808730f85c39ee2cbeaed2c6217e9d546cd445b3bc897f"
OCTOBER_SOURCE_MANIFEST_SHA256 = "b9a38083131dcf1a4c7c4e73ebe9a97dff2a4b7ef91a8df53f8bc38cadfa48db"
OCTOBER_LOCAL_MANIFEST_SHA256 = "fd46f09e642fab95dbd57ddc8a3f79461ff31c9bb9c21772b7e4729effee1fd2"
OCTOBER_RECEIPT_COMMIT = "2cc0aa74f5d05d1a7b500899225d0c1c148d957f"
OCTOBER_SPATIAL_RECEIPT_SHA256 = "743abf565518d588fc434d1cff74f58a507f46c7eee4b25df92222f21b330ba6"
OCTOBER_FINAL_RECEIPT_URL = f"https://github.com/pourpoint-tools/hfx/blob/{OCTOBER_RECEIPT_COMMIT}/docs/v0.4.0/receipts/grit-final-delivery.json"
OCTOBER_FINAL_RECEIPT_SHA256 = "374f4f100106e600b608a02fa8f3cd0939023ec73ec72c69ee425d9149a08f26"
OCTOBER_ORIGINAL_VERIFICATION_URL = (
    f"https://github.com/pourpoint-tools/hfx/blob/{OCTOBER_RECEIPT_COMMIT}/docs/v0.4.0/grit-river-delivery.md"
)
OCTOBER_ORIGINAL_VERIFICATION_SHA256 = "0bb03e219657775da72271573105dd8d2451a4d9def9c2c1cd8b50d00a4be9aa"
PRESERVED_OCTOBER_MANIFEST_PATH = f"provenance/manifest-{OCTOBER_MANIFEST_SHA256}.json"
CONNECTIVITY_PROVENANCE_PATH = "provenance/connectivity-completeness-2026-11-23.json"

EXPECTED_PREVIOUS_CONNECTIVITY: list[dict[str, object]] = [
    {"level": 0, "status": "not_supplied"},
    {"level": 1, "status": "supplied"},
    {"level": 2, "status": "supplied"},
]

UPDATED_CONNECTIVITY: list[dict[str, object]] = [
    {"level": 0, "status": "not_supplied"},
    {"level": 1, "status": "supplied", "completeness": "complete"},
    {"level": 2, "status": "supplied", "completeness": "complete"},
]

_MODIFIED_TOP_LEVEL_KEYS = frozenset(
    {"catchment_connectivity", "relationships", "artifacts", "connectivity_completeness_provenance"}
)


def fingerprint_bytes(data: bytes) -> Fingerprint:
    """Compute size, full SHA-256, and bounded head/tail SHA-256 digests."""
    size = len(data)
    return {
        "size_bytes": size,
        "sha256": hashlib.sha256(data).hexdigest(),
        "head_sha256": hashlib.sha256(data[:65536]).hexdigest(),
        "tail_sha256": hashlib.sha256(data[max(0, size - 65536) :]).hexdigest(),
    }


@dataclass(frozen=True)
class PreparedMetadataUpdate:
    """In-memory staged artifacts for updating GRIT connectivity completeness metadata."""

    preserved_manifest_path: str
    preserved_manifest_bytes: bytes
    provenance_path: str
    provenance_bytes: bytes
    updated_manifest_bytes: bytes


def _parse_manifest_object(raw: bytes, label: str) -> dict[str, object]:
    loaded = json.loads(raw.decode("utf-8"))
    if not isinstance(loaded, dict):
        raise ValueError(f"{label} must decode to a JSON object")
    return cast(dict[str, object], loaded)


def _validate_base_manifest(manifest: dict[str, object]) -> None:
    if manifest.get("format_version") != "0.4.0":
        raise ValueError("expected HFX format_version 0.4.0")
    if manifest.get("representation") != "spatial-geoparquet-adjacency" or manifest.get("representation_version") != 1:
        raise ValueError("expected spatial-geoparquet-adjacency representation version 1")
    if manifest.get("fabric_name") != "grit" or manifest.get("fabric_version") != "1.0":
        raise ValueError("expected GRIT 1.0 manifest")
    if manifest.get("crs") != "EPSG:4326":
        raise ValueError("expected EPSG:4326 CRS")
    if manifest.get("input_manifest_sha256") != OCTOBER_SOURCE_MANIFEST_SHA256:
        raise ValueError("unexpected input_manifest_sha256")


def _build_provenance_payload(
    previous: dict[str, object],
    *,
    previous_manifest_bytes: bytes,
    previous_sha256: str,
    preserved_manifest_path: str,
) -> dict[str, object]:
    families = cast(dict[str, dict[str, object]], previous["families"])
    artifacts = cast(dict[str, object], previous["artifacts"])
    external_artifacts = cast(dict[str, object], previous.get("external_artifacts", {}))
    return {
        "bundle_url": PUBLIC_BUNDLE_URL,
        "bundle_identity": {
            "format_version": previous["format_version"],
            "representation": previous["representation"],
            "representation_version": previous["representation_version"],
            "fabric_name": previous["fabric_name"],
            "fabric_version": previous["fabric_version"],
            "region": previous["region"],
            "unit_count": previous["unit_count"],
            "rivers_count": families["rivers"]["count"],
            "river_nodes_count": families["river_nodes"]["count"],
            "input_manifest_sha256": previous["input_manifest_sha256"],
            "source_manifest": previous["source_manifest"],
            "vector_prefix": "vectors/2026-10-18-final/",
        },
        "preserved_previous_manifest": {
            "path": preserved_manifest_path,
            "url": PUBLIC_BUNDLE_URL + preserved_manifest_path,
            "sha256": previous_sha256,
            "size_bytes": len(previous_manifest_bytes),
        },
        "october_verification": {
            "commit": OCTOBER_RECEIPT_COMMIT,
            "original_verification_url": OCTOBER_ORIGINAL_VERIFICATION_URL,
            "original_verification_sha256": OCTOBER_ORIGINAL_VERIFICATION_SHA256,
            "final_receipt_url": OCTOBER_FINAL_RECEIPT_URL,
            "final_receipt_sha256": OCTOBER_FINAL_RECEIPT_SHA256,
            "public_manifest_sha256": OCTOBER_MANIFEST_SHA256,
            "source_manifest_sha256": OCTOBER_SOURCE_MANIFEST_SHA256,
            "local_manifest_sha256": OCTOBER_LOCAL_MANIFEST_SHA256,
            "local_spatial_verification_receipt": {
                "artifact": "final-spatial-verification.json",
                "sha256": OCTOBER_SPATIAL_RECEIPT_SHA256,
                "size_bytes": 2134,
            },
            "verified_counts": {
                "catchments": 22419857,
                "rivers": 22337300,
                "river_nodes": 22278858,
                "catchment_flow_edges": 22395952,
                "river_flow_edges": 22395952,
                "local_river_catchment_pairs": 22337300,
            },
            "summary": (
                "Original GRIT verifier paired every segment and reach polygon with its same-tier "
                "source river and checked delivered catchment edges against source edges in both "
                "directions; final spatial conversion verified all immediate relationships against "
                "that delivery."
            ),
        },
        "catchment_connectivity": [
            {
                "level": 0,
                "status": "not_supplied",
                "source_meaning": "GRIT supplies component membership only at this level.",
            },
            {
                "level": 1,
                "status": "supplied",
                "completeness": "complete",
                "source_meaning": "All immediate connections from the GRIT v1.0 segment network are represented.",
            },
            {
                "level": 2,
                "status": "supplied",
                "completeness": "complete",
                "source_meaning": "All immediate connections from the GRIT v1.0 reach network are represented.",
            },
        ],
        "reused_artifacts": {
            "vector_artifact_count": len(artifacts),
            "external_raster_keys": sorted(external_artifacts),
            "all_vector_and_raster_objects_unchanged": True,
        },
        "completeness_scope": (
            "Completeness is relative to GRIT v1.0's published vector network (Zenodo record 17435232, "
            "CC BY-NC 4.0). It does not certify vector-raster consistency "
            "(https://github.com/pourpoint-tools/hfx/issues/264); known GRIT v1.0 drainage-direction "
            "raster merge loops are documented separately "
            "(https://github.com/pourpoint-tools/hfx/issues/269#issuecomment-7104717277). "
            "The existing HFX 0.4.0 contract supports complete levels without a per-catchment cut column."
        ),
        "references": {
            "hfx_issue": "https://github.com/pourpoint-tools/hfx/issues/288",
            "hfx_vision": "planning/visions/2026-11-23-grit-catchment-connectivity-completeness.md",
            "pourpoint_issue": "https://github.com/pourpoint-tools/pourpoint/issues/227",
            "vector_raster_consistency_issue": "https://github.com/pourpoint-tools/hfx/issues/264",
            "grit_raster_loop_finding": "https://github.com/pourpoint-tools/hfx/issues/269#issuecomment-7104717277",
        },
    }


def prepare_metadata_update(
    previous_manifest_bytes: bytes,
    *,
    expected_previous_sha256: str = OCTOBER_MANIFEST_SHA256,
    provenance_path: str = CONNECTIVITY_PROVENANCE_PATH,
) -> PreparedMetadataUpdate:
    """Validate the October GRIT manifest and stage the connectivity-completeness metadata update."""
    actual_sha256 = hashlib.sha256(previous_manifest_bytes).hexdigest()
    if actual_sha256 != expected_previous_sha256:
        raise ValueError(f"previous manifest digest {actual_sha256} does not match expected {expected_previous_sha256}")
    previous = _parse_manifest_object(previous_manifest_bytes, "previous manifest")
    _validate_base_manifest(previous)
    if previous.get("catchment_connectivity") != EXPECTED_PREVIOUS_CONNECTIVITY:
        raise ValueError("unexpected previous catchment_connectivity declaration")
    relationships = cast(dict[str, dict[str, object]], previous.get("relationships", {}))
    catchment_rel = relationships.get("catchments")
    if not isinstance(catchment_rel, dict) or catchment_rel.get("connectivity") != EXPECTED_PREVIOUS_CONNECTIVITY:
        raise ValueError("unexpected previous catchment relationship connectivity declaration")
    if catchment_rel.get("connectivity_source") != "manifest.json":
        raise ValueError("unexpected previous catchment relationship connectivity_source")

    preserved_manifest_path = f"provenance/manifest-{expected_previous_sha256}.json"
    previous_artifacts = cast(dict[str, object], previous.get("artifacts", {}))
    for path in (preserved_manifest_path, provenance_path):
        if path in previous_artifacts:
            raise ValueError(f"artifact path already exists in previous manifest: {path}")

    provenance_payload = _build_provenance_payload(
        previous,
        previous_manifest_bytes=previous_manifest_bytes,
        previous_sha256=expected_previous_sha256,
        preserved_manifest_path=preserved_manifest_path,
    )
    provenance_bytes = (json.dumps(provenance_payload, indent=2) + "\n").encode("utf-8")

    updated = copy.deepcopy(previous)
    updated["catchment_connectivity"] = copy.deepcopy(UPDATED_CONNECTIVITY)
    cast(dict[str, dict[str, object]], updated["relationships"])["catchments"]["connectivity"] = copy.deepcopy(
        UPDATED_CONNECTIVITY
    )
    updated["connectivity_completeness_provenance"] = {
        "artifact": provenance_path,
        "previous_manifest": preserved_manifest_path,
        "previous_manifest_sha256": expected_previous_sha256,
        "october_final_receipt": OCTOBER_FINAL_RECEIPT_URL,
        "october_original_verification": OCTOBER_ORIGINAL_VERIFICATION_URL,
        "scope": (
            "Relative to GRIT v1.0 published vector network (Zenodo 17435232, CC BY-NC 4.0); "
            "does not certify vector-raster consistency (HFX #264) or alter known GRIT v1.0 "
            "drainage-direction raster loops (HFX #269)."
        ),
    }
    combined_artifacts = dict(previous_artifacts)
    combined_artifacts[preserved_manifest_path] = fingerprint_bytes(previous_manifest_bytes)
    combined_artifacts[provenance_path] = fingerprint_bytes(provenance_bytes)
    updated["artifacts"] = {key: combined_artifacts[key] for key in sorted(combined_artifacts)}
    updated_manifest_bytes = (json.dumps(updated, indent=2) + "\n").encode("utf-8")

    return PreparedMetadataUpdate(
        preserved_manifest_path=preserved_manifest_path,
        preserved_manifest_bytes=previous_manifest_bytes,
        provenance_path=provenance_path,
        provenance_bytes=provenance_bytes,
        updated_manifest_bytes=updated_manifest_bytes,
    )


def verify_metadata_update(
    previous_manifest_bytes: bytes,
    updated_manifest_bytes: bytes,
    preserved_manifest_bytes: bytes,
    provenance_bytes: bytes,
    *,
    expected_previous_sha256: str = OCTOBER_MANIFEST_SHA256,
    provenance_path: str = CONNECTIVITY_PROVENANCE_PATH,
) -> dict[str, object]:
    """Independently verify that only connectivity completeness and provenance were updated."""
    if hashlib.sha256(previous_manifest_bytes).hexdigest() != expected_previous_sha256:
        raise ValueError("previous manifest digest mismatch")
    if preserved_manifest_bytes != previous_manifest_bytes:
        raise ValueError("preserved manifest bytes differ from previous manifest")
    if hashlib.sha256(preserved_manifest_bytes).hexdigest() != expected_previous_sha256:
        raise ValueError("preserved manifest digest mismatch")

    previous = _parse_manifest_object(previous_manifest_bytes, "previous manifest")
    updated = _parse_manifest_object(updated_manifest_bytes, "updated manifest")
    _validate_base_manifest(updated)

    preserved_manifest_path = f"provenance/manifest-{expected_previous_sha256}.json"

    prev_other = {k: v for k, v in previous.items() if k not in _MODIFIED_TOP_LEVEL_KEYS}
    upd_other = {k: v for k, v in updated.items() if k not in _MODIFIED_TOP_LEVEL_KEYS}
    if prev_other != upd_other:
        raise ValueError("unexpected top-level manifest changes outside connectivity and provenance")

    if updated.get("catchment_connectivity") != UPDATED_CONNECTIVITY:
        raise ValueError("updated catchment_connectivity does not match expected declarations")
    connectivity_entries = cast(list[dict[str, object]], updated["catchment_connectivity"])
    for entry in connectivity_entries:
        if entry["status"] == "not_supplied":
            if any(key in entry for key in ("completeness", "upstream_cuts_column", "derivation")):
                raise ValueError("not_supplied catchment_connectivity entry must omit completeness")
        elif entry.get("completeness") != "complete" or "upstream_cuts_column" in entry:
            raise ValueError("supplied catchment_connectivity entry must declare complete without cut column")

    prev_rels = cast(dict[str, dict[str, object]], previous["relationships"])
    upd_rels = cast(dict[str, dict[str, object]], updated["relationships"])
    if set(prev_rels) != set(upd_rels):
        raise ValueError("relationship domains changed")
    for domain in prev_rels:
        if domain != "catchments" and prev_rels[domain] != upd_rels[domain]:
            raise ValueError(f"non-catchment relationship metadata changed for {domain}")

    prev_cat_rel = dict(prev_rels["catchments"])
    upd_cat_rel = dict(upd_rels["catchments"])
    if upd_cat_rel.get("connectivity") != updated["catchment_connectivity"]:
        raise ValueError("catchment relationship connectivity does not match root catchment_connectivity")
    if upd_cat_rel.get("connectivity_source") != "manifest.json":
        raise ValueError("catchment relationship connectivity_source changed")
    prev_cat_rel.pop("connectivity", None)
    upd_cat_rel.pop("connectivity", None)
    if prev_cat_rel != upd_cat_rel:
        raise ValueError("catchment relationship non-connectivity fields changed")

    prev_artifacts = cast(dict[str, object], previous["artifacts"])
    upd_artifacts = cast(dict[str, object], updated["artifacts"])
    expected_added = {preserved_manifest_path, provenance_path}
    if set(upd_artifacts) != set(prev_artifacts) | expected_added:
        raise ValueError("updated artifacts keys differ from previous plus new provenance artifacts")
    for key, fp in prev_artifacts.items():
        if upd_artifacts[key] != fp:
            raise ValueError(f"existing artifact fingerprint must remain unchanged: {key}")
    if upd_artifacts[preserved_manifest_path] != fingerprint_bytes(preserved_manifest_bytes):
        raise ValueError("preserved manifest fingerprint mismatch")
    if upd_artifacts[provenance_path] != fingerprint_bytes(provenance_bytes):
        raise ValueError("provenance artifact fingerprint mismatch")

    prov_decl = updated.get("connectivity_completeness_provenance")
    if not isinstance(prov_decl, dict):
        raise ValueError("missing connectivity_completeness_provenance declaration")
    if (
        prov_decl.get("artifact") != provenance_path
        or prov_decl.get("previous_manifest") != preserved_manifest_path
        or prov_decl.get("previous_manifest_sha256") != expected_previous_sha256
    ):
        raise ValueError("connectivity_completeness_provenance declaration references mismatch")

    provenance = _parse_manifest_object(provenance_bytes, "provenance artifact")
    october = provenance.get("october_verification")
    preserved_meta = provenance.get("preserved_previous_manifest")
    if not isinstance(october, dict) or not isinstance(preserved_meta, dict):
        raise ValueError("provenance artifact missing October verification or preserved manifest records")
    if october.get("final_receipt_url") != OCTOBER_FINAL_RECEIPT_URL:
        raise ValueError("provenance artifact missing October final receipt URL")
    if october.get("original_verification_url") != OCTOBER_ORIGINAL_VERIFICATION_URL:
        raise ValueError("provenance artifact missing October original verification URL")
    if (
        preserved_meta.get("sha256") != expected_previous_sha256
        or preserved_meta.get("path") != preserved_manifest_path
    ):
        raise ValueError("provenance artifact preserved manifest metadata mismatch")
    if not isinstance(provenance.get("completeness_scope"), str) or not provenance["completeness_scope"]:
        raise ValueError("provenance artifact missing completeness_scope")

    return {
        "status": "passed",
        "previous_manifest_sha256": expected_previous_sha256,
        "updated_manifest_sha256": hashlib.sha256(updated_manifest_bytes).hexdigest(),
        "updated_manifest_size_bytes": len(updated_manifest_bytes),
        "preserved_manifest_path": preserved_manifest_path,
        "provenance_path": provenance_path,
        "provenance_sha256": hashlib.sha256(provenance_bytes).hexdigest(),
        "reused_artifact_count": len(prev_artifacts),
    }


def _write_candidate_directory(output_dir: Path, prepared: PreparedMetadataUpdate) -> None:
    if output_dir.exists():
        raise FileExistsError(output_dir)
    output_dir.mkdir(parents=True)
    for rel_path, content in (
        ("manifest.json", prepared.updated_manifest_bytes),
        (prepared.preserved_manifest_path, prepared.preserved_manifest_bytes),
        (prepared.provenance_path, prepared.provenance_bytes),
    ):
        dest = output_dir / rel_path
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(content)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare", help="Stage updated manifest and provenance artifacts")
    prepare_parser.add_argument("--input-manifest", type=Path, required=True)
    prepare_parser.add_argument("--output-dir", type=Path, required=True)
    prepare_parser.add_argument("--expected-previous-sha256", default=OCTOBER_MANIFEST_SHA256)
    prepare_parser.add_argument("--provenance-path", default=CONNECTIVITY_PROVENANCE_PATH)

    verify_parser = subparsers.add_parser("verify", help="Verify staged manifest and provenance artifacts")
    verify_parser.add_argument("--previous-manifest", type=Path, required=True)
    verify_parser.add_argument("--candidate-dir", type=Path, required=True)
    verify_parser.add_argument("--expected-previous-sha256", default=OCTOBER_MANIFEST_SHA256)
    verify_parser.add_argument("--provenance-path", default=CONNECTIVITY_PROVENANCE_PATH)

    args = parser.parse_args()
    if args.command == "prepare":
        previous_bytes = args.input_manifest.read_bytes()
        prepared = prepare_metadata_update(
            previous_bytes,
            expected_previous_sha256=args.expected_previous_sha256,
            provenance_path=args.provenance_path,
        )
        summary = verify_metadata_update(
            previous_bytes,
            prepared.updated_manifest_bytes,
            prepared.preserved_manifest_bytes,
            prepared.provenance_bytes,
            expected_previous_sha256=args.expected_previous_sha256,
            provenance_path=args.provenance_path,
        )
        _write_candidate_directory(args.output_dir, prepared)
        print(json.dumps(summary, indent=2))
        return

    previous_bytes = args.previous_manifest.read_bytes()
    preserved_path = f"provenance/manifest-{args.expected_previous_sha256}.json"
    summary = verify_metadata_update(
        previous_bytes,
        (args.candidate_dir / "manifest.json").read_bytes(),
        (args.candidate_dir / preserved_path).read_bytes(),
        (args.candidate_dir / args.provenance_path).read_bytes(),
        expected_previous_sha256=args.expected_previous_sha256,
        provenance_path=args.provenance_path,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
