# GRIT catchment connectivity completeness

Ship the public GRIT HFX 0.4.0 bundle alongside Pourpoint 0.5.0 at the URL
already used by Pourpoint’s Rhine-at-Basel example. A first-time user should be
able to select the upstream segment catchments with default settings and receive
an explicit `complete` coverage result. Today the bundle supplies their flow
graph but omits the completeness declaration, so Pourpoint requires opt-in
even for this ordinary selection. [Public manifest](https://basin-delineations-public.upstream.tech/grit/hfx-v0.4.0/manifest.json) · [Pourpoint example](https://github.com/pourpoint-tools/pourpoint/blob/28b4ea153dd3ab7dfb300de3827e61dc1951c766/README.md#delineate-the-rhine-at-basel)

## What the October verification establishes

Use the October worldwide run as the evidence for this metadata-only update.
The original verifier paired each segment and reach polygon with its same-tier
source river and checked the delivered catchment edges against the source edges
in both directions. The final spatial conversion verified the immediate
relationships against that delivery. Its receipt identifies the public manifest
by SHA-256 `04bbfb50907feaf339808730f85c39ee2cbeaed2c6217e9d546cd445b3bc897f`.
[Original verification](https://github.com/pourpoint-tools/hfx/blob/2cc0aa74f5d05d1a7b500899225d0c1c148d957f/docs/v0.4.0/grit-river-delivery.md) · [Final receipt](https://github.com/pourpoint-tools/hfx/blob/2cc0aa74f5d05d1a7b500899225d0c1c148d957f/docs/v0.4.0/receipts/grit-final-delivery.json)

| Catchment level | Connectivity declaration | Source meaning |
| --- | --- | --- |
| 0 · component | `status: not_supplied`, omit completeness | GRIT supplies component membership only at this level. |
| 1 · segment | `status: supplied`, `completeness: complete` | All immediate connections from the GRIT v1.0 segment network are represented. |
| 2 · reach | `status: supplied`, `completeness: complete` | All immediate connections from the GRIT v1.0 reach network are represented. |

Apply these declarations consistently in the root manifest and its catchment
relationship metadata. Record the October receipt, bundle identity, and scope
of the completeness claim in bundle provenance. Preserve the October manifest
by hash and reuse the verified vector and raster objects. The existing HFX
0.4.0 contract supports complete levels without a per-catchment cut column.
[Bundle organization](../../docs/site/bundle.md) · [Vector schemas](../../docs/site/vector-data.md)

Completeness is relative to GRIT’s published vector network. Keep the known
GRIT v1.0 drainage-direction loops documented separately. Verifying agreement
between a polygon watershed and raster drainage belongs to the vector–raster
consistency effort. [GRIT raster finding](https://github.com/pourpoint-tools/hfx/issues/269#issuecomment-7104717277) · [HFX #264](https://github.com/pourpoint-tools/hfx/issues/264)

## Release verification

Verify a candidate with the generic coverage checks already delivered in
[Pourpoint #227](https://github.com/pourpoint-tools/pourpoint/issues/227).
Use the release’s normal selection API and default coverage policy.

- Follow the Rhine-at-Basel example through snapping and upstream selection at
  its default segment level. Confirm the selection contains the expected upstream
  catchments and reports complete coverage.
- Select upstream catchments at reach level for a known connected GRIT outlet.
  Confirm the selection includes the expected upstream reaches’ catchments and
  reports complete coverage.
- Request upstream catchment selection at component level. Confirm Pourpoint
  reports `selection.connectivity_unavailable`.

Coordinate publication of the verified metadata at the existing public URL in
time for Pourpoint 0.5.0. Repeat all three checks against that URL with a fresh
cache and record the public manifest digest, Pourpoint revision, targets, and
observed results. Completing this delivery closes [HFX #288](https://github.com/pourpoint-tools/hfx/issues/288).
Publishing this vision leaves the issue open. Paid infrastructure requires
separate authorization.
