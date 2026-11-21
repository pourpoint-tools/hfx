# Vector feature schemas

These schemas apply to the [current representation](bundle.md). Physical access
is specified in [spatial layout](spatial-layout.md).

## Catchments

The catchment GeoParquet file contains one row per drainage unit.
Attributes and WKB occupy the same row. The required fields are:

| Field | Arrow type | Nullable | Meaning |
| --- | --- | --- | --- |
| `id` | `int64` | No | Positive dataset-local unique unit ID |
| `level` | `int16` | No | Nonnegative dataset-local resolution tier |
| `parent_id` | `int64` | Yes | Containing unit at a strictly coarser available level |
| `area_km2` | `float32` | No | Local geodesic area in km² |
| `up_area_km2` | `float32` | Yes | Inclusive same-level upstream area in km², under documented dataset semantics |
| `outlet_lon` | `float64` | Yes | Finite EPSG:4326 longitude of a single source-described outlet |
| `outlet_lat` | `float64` | Yes | Finite EPSG:4326 latitude of that outlet |
| `bbox` | `struct<xmin: float32!, ymin: float32!, xmax: float32!, ymax: float32!>` | No | Finite ordered, outward-covering EPSG:4326 bounds; `!` denotes non-nullable leaves |
| `geometry` | `binary` | No | Exact Polygon or MultiPolygon WKB in EPSG:4326 |

Source identity is retained as `source_namespace: string` and `source_id: string`
(non-nullable). A namespace identifies the source product and layer, including
region/basin where needed for unique source identity. The pair identifies a source
feature; preserve native ID columns and types independently of the string rendering.
Only unchanged catchment artifacts from the earlier normalized 0.4.0 publication
may retain their existing optional source identity fields when documented provenance
establishes the mapping. Newly produced catchments require both nonnullable fields.

Smaller levels are coarser. Each coarsest unit has null `parent_id`; other nonnull
parents resolve at a strictly coarser available catchment level. A child polygon
must lie within its declared parent, and sibling polygon interiors must not overlap,
within explicitly documented numerical geometry tolerances. A producer declares
the tolerance values, units and comparison method and reports actual validation
coverage. This contract does not require exhaustive worldwide polygon topology
certification; full reference and source-fidelity verification is a separate duty.
Tolerances account for numerical representation, not arbitrary source discrepancies. A violation beyond
them is nonconforming; recording an anomaly does not waive the invariant. Investigate
such originals and obtain an explicit contract decision before production rather
than silently repairing geometry or relaxing containment. Record validation
coverage separately: sampled checks do not establish worldwide conformance.

For MERIT-Basins source catchments only, an explicitly approved narrow ring-normalization
rule admits invalid source polygons whose raw OGC validity failure begins with
`Ring Self-intersection[` (raster-to-polygon self-touching rings at pixel corners),
without authorizing arbitrary geometry cleanup, simplification, hole filling, or a new
HFX format version. Already valid source geometries remain unchanged. An eligible
MERIT repair is admitted only when every check below passes without relaxing
thresholds:

1. Output is a valid, non-empty `Polygon` or `MultiPolygon` without Z or M ordinates
   and with finite WGS84 coordinates within declared bounds.
2. Exact IEEE 754 `float64` source vertex set is preserved.
3. Exact IEEE 754 `float64` source undirected edge multiset is preserved, retaining all
   holes and disconnected components without moving, filling, simplifying, or
   discarding source land.
4. WGS84 geodesic area changes by at most `max(1e-9, 1e-12 * raw_area_km2)` km².
5. Sibling interior non-overlap (`DE-9IM` pattern `2********` is false) holds against all
   relevant delivered neighbors globally across partitions and coastal hillslopes,
   including other candidate repairs; pre-existing overlap is not waived.

Features outside this rule or failing any acceptance check remain omitted with explicit
reasons in `omitted-features.json`. Real-area MERIT catchment gaps are restored only
through valid source-backed polygons under this rule and are never bridged by graph
contraction, preserving TDX-Hydro's distinct zero-length, zero-area virtual-link
contraction rule without broadening it. Every MERIT catchment records non-nullable
`geometry_normalization: string` (`"original"` or `"ring_self_intersection_split"`,
declared in `categorical_columns`), per-feature audit records with original WKB, raw
failure reason, ring/part counts, algorithm, and measured checks in
`repaired-features.jsonl.gz`, and library versions and exhaustive per-partition counts
in `source-metadata.json`.

Components may have multiple or no represented outlets, including coastal and
inland sink polygon-only units. Outlet coordinates are either both null or both
present. Supply them only for a single source-described outlet; never choose an
arbitrary representative from multiple outlets or invent a point. Source river
nodes and associations retain multiple outlet facts independently.

`catchment_connectivity` distinguishes supplied edges from absent information.
For a `not_supplied` level, every graph row has an empty upstream list; that is not
an assertion of natural termination. The topology declaration applies only to
supplied levels. Even there, missing links alone do not classify a terminal.
Preserve source explanations and distinguish subset cuts. Additional typed source
attributes and nulls remain intact. Preserve native local polygon area separately
from geodesic `area_km2`; partitioned drainage area, mainstem-following area and
upstream polygon-union area are different quantities. Do not relabel one as another.
Document upstream-area semantics, derivations and anomalies; use null rather than
an unsupported value.

For supplied catchment levels with `completeness: "partial"`, producers may supply
per-unit completeness and boundary fields:

- `upstream_boundary_cuts` (`int32`, nullable): nonnegative count of immediate
  incoming source edges whose upstream catchment is absent or omitted (`0` when all
  immediate incoming source edges resolve to delivered catchments; `> 0` when
  immediate incoming edges were cut). `null` means incoming catchment connectivity
  could not be verified from valid source topology (such as unchannelized coastal
  hillslopes without a topology table or catchments whose paired source river was
  omitted for degenerate geometry).
- Optional dataset `downstream_boundary_status` (`string`, non-nullable, documented in
  `categorical_columns`) distinguishing `"internal"` downstream edges from source zero
  sentinels, downstream cuts to missing catchments, omitted degenerate rivers, and
  unassociated coastal hillslopes.

## Rivers

| Field | Arrow type | Nullable | Meaning |
| --- | --- | --- | --- |
| `id` | `int64` | No | Positive dataset-local river ID, unique across river levels |
| `level` | `int16` | No | Nonnegative river tier; smaller is coarser |
| `parent_id` | `int64` | Yes | Constituent membership in a river at a strictly coarser available tier |
| `source_namespace`, `source_id` | `string` | No | Source layer identity and ID rendering; retain native ID separately |
| `local_catchment_id` (optional column) | `int64` | Yes | Sole represented local drainage association; see plural association rules |
| `component_catchment_id` | `int64` | Yes | Associated drainage unit representing the connected network component |
| `upstream_node_id`, `downstream_node_id` | `int64` | Yes | Source-directed endpoint references to node rows |
| `geometry_orientation` | `int8` | No | 0 unspecified, 1 upstream-first, 2 downstream-first |
| `bbox` | Catchment bbox struct | No | Outward-covering EPSG:4326 bounds |
| `geometry` | `binary` | No | LineString or MultiLineString WKB in EPSG:4326 |

The `local_catchment_id` column is optional. The authoritative plural
[local drainage associations](spatial-layout.md#local-drainage-associations)
determine whether it may be included and what its values mean.

Within each feature table, `(source_namespace, source_id)` is unique across all
levels; river, node and catchment tables have separate HFX ID domains. Namespaces
retain source layer identity, not a label inferred from coincident geometry.
River IDs are not row numbers or globally stable IDs across datasets/versions.
River and catchment tiers are independent. One-tier rivers need no invented
parents. Stream order is an attribute, not a tier or universal snapping weight.
Parent membership does not require equality with the union of child geometries.
Each nonnull catchment association resolves to a catchment ID; multiple rivers may
share either target. Sources need not supply components or polygons. Null means no
represented association, with absence versus unavailable evidence documented.
Never infer hydrological association solely from spatial intersection.

Geometry retains coordinate and part order. Do not silently reverse, simplify,
clip or repair it. Record necessary transformations, including source geometry
export and reprojection. Orientation 1/2 asserts first-to-last order relative to
flow; for MultiLineString it also asserts parts are in connected flow order with
each part consistently oriented. Otherwise use 0. Connectivity and directed node
references remain authoritative regardless of orientation; 0 is not a failure.

## Source nodes

The optional `river_nodes` family has no flow graph or topology declaration. Required node fields are nonnull `id: int64`, `level: int16`,
`source_namespace: string`, `source_id: string`, `bbox` of the catchment struct
type, and `geometry: binary` (Point WKB in EPSG:4326), plus nullable
`coarser_node_id: int64`. Node IDs are positive and unique across node tiers;
`level` is the nonnegative river tier. A coarser-node reference resolves at a strictly smaller
level and means source-declared correspondence, not constituent membership.
Keep source layers distinct even if coordinates and numeric source IDs match.

Nonnull river endpoint IDs resolve to nodes at that river's level. These references
are the canonical incidence relation: incoming rivers have that downstream node;
outgoing rivers have that upstream node. Preserve original incidence strings,
ordering-dependent width shares, cross-tier reference columns, diagnostics and
other typed attributes alongside it. Do not invent nodes for sources without them.
No node reference may be nonnull when the node artifact is absent. Verify that
source-supplied flow edges and endpoint incidence agree; sharing a node alone does
not imply every possible incoming/outgoing river pair is a source flow edge.

## Flow connectivity

Supplied flow graphs are acyclic. `tree` allows at most one downstream neighbor;
`dag` allows splits and reconvergence without duplicating features. Each upstream
ID denotes an edge from that feature to the row's feature at the same level.
IDs resolve within the corresponding domain. Lists have no nulls, duplicates or
self-links. Parent membership is not flow. Validate complete source edges before
declaring conformance; a cycle cannot be excused by an anomaly note or removed
silently. Acyclicity does not certify hydrological continuity.

When catchment connectivity is derived from supplied river connectivity and verified
1-to-1 local drainage associations, only immediate edges whose upstream and
downstream rivers both have delivered associated catchments are materialized in the
catchment flow graph; omitted or missing polygons with drainage area (such as TDX
`LINKNO == 0` omissions) are never contracted or bridged across gaps, and missing
links are never turned into natural outlets.

Where a river network represents multi-way confluences with zero-length, zero-area
virtual links (for example, TauDEM binary junctions in TDX-Hydro where `Length == 0`,
`geom.length == 0`, `DSContArea == USContArea`, `upstream_degree >= 1`, and
`downstream_degree == 1`, failing closed on any null or non-finite value),
exhaustively verified virtual stubs may be contracted in the derived catchment flow
graph while leaving the river flow graph untouched. Transitive contraction propagates
both upstream catchments and incoming boundary cuts into the downstream catchment and
never bridges across a non-contracted polygonless river.

When a level declares `completeness: "partial"` with `upstream_cuts_column`, an
upstream catchment selection is complete relative to the source network if and only
if every visited unit in the traversed closure has non-null
`upstream_boundary_cuts == 0` (or the declared column equals `0`). If any visited
unit has a positive cut count or `null`, consumers must not return an ordinary
complete result; partial results require explicit opt-in and missing-coverage
reporting. Provenance metadata alone does not enforce this runtime check; consumer
support requires explicit coordination.

## Typed meanings and source fidelity

HFX fixes structural fields, not a complete source attribute schema. Preserve
useful source attributes, physical types, nulls and meanings. Explicitly approved
lossless type mappings may change physical width while preserving exact values and
nulls; retain original per-source declarations and document the mapping. TDX
catchment `streamID` widens int32 to int64 under this rule. Where a source exports
multi-part catchments as separate single-part `POLYGON` rows sharing an association
key (TDX basins `1020018110` and `5020049720`), collecting those fragments into one
OGC-valid `MultiPolygon` without dissolving or altering coordinates retains
`min(source_fids)` in `fid` and `source_id` (keeping `source_id` unique within
`source_namespace`) and preserves every original `fid` in ascending order in nonnull
`source_fids: list<int64>`. Document retained,
renamed, transformed and omitted fields; a collision with an HFX structural name
requires an explicit source-column rename map in provenance, not value replacement.
Retain source field/schema metadata where applicable; update geometry metadata for
actual delivered CRS and bounds rather than carrying stale projected metadata.

The optional `categorical_columns` array binds a dictionary to one scalar column:

```json
{
  "artifact": "rivers.parquet", "column": "level", "physical_type": "int16",
  "codes": [{"code": 0, "meaning": "segment"}, {"code": 1, "meaning": "reach"}],
  "null_meaning": "not permitted"
}
```

`artifact` is an exact local artifact path, `column` a top-level field name.
`physical_type` is one of `string`, `bool`, `int8`, `int16`, `int32`, `int64`,
`uint8`, `uint16`, `uint32`, `uint64`, matching the Arrow field's logical scalar
type (dictionary encoding, if present, does not change this type). Codes are unique,
nonnull JSON values of that type; integer parsing must preserve full precision.
Every observed nonnull value has one nonempty string meaning. `null_meaning` is a
required nonempty description of null interpretation, even for nonnullable fields.
At most one dictionary binds each artifact/column pair. Dictionaries define dataset
semantics, not cross-dataset equivalence or Arrow compression. Level descriptions
supplement HFX ordering, never replace parents/edges. A source `termination_kind`
may use this mechanism; no universal outlet vocabulary is required. Preserve
underspecified source codes with an explicit uncertainty description.

Subsets retain external source edge identities in source columns and document
which edges were cut, in each direction; structural graph references still resolve
within the delivered subset. Never turn a known cut into a source terminal.

