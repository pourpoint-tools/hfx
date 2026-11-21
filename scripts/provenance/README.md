# Retained production provenance

The GRIT producer and its independent verification tests consume
`grit-original-inventory.json` and `grit-raster-delivery.json` as pinned source
identities and raster declarations. These records are retained production inputs,
not reader-facing documentation or new validation results. Do not replace their
identities with hashes from an untrusted download.

`NOTICE` and `CITATION.txt` preserve the source attribution formerly carried by
the raster experiment. GRIT data is CC BY-NC 4.0; repository code is MIT.
Historical build, delivery and experiment records remain in Git history.

`merit-delivery.json`, `merit-NOTICE`, and `merit-CITATION.txt` record the pinned
source inventories/hashes, global normalization receipts, private Hetzner S3
artifact identities, remote range-probe results, spend accounting, and license
attributions for the MERIT HFX 0.4.0 delivery (`#281`). They are secret-free
provenance records for `MERIT Hydro` (CC BY-NC 4.0) and `MERIT-Basins`
(CC BY-NC-SA 4.0 / UCSD non-commercial research terms).

`merit-catchment-connectivity.json` records the secret-free provenance for the
candidate MERIT HFX 0.4.0 bundle (`#285`, disjoint candidate prefix) derived
from the verified `#281` base bundle with derived `catchment_connectivity`,
mixed raw-shapefile/base-bundle verification receipts, remote range-probe and
isolated consumer verification results, and bounded Hetzner Cloud spend.

`merit-geometry-restoration.json` records the secret-free provenance for the
promoted MERIT HFX 0.4.0 bundle
(`s3://pourpoint-hfx/merit/hfx-0.4.0/merit-hydro-v1.0.1-basins-v1.0-bugfix1-20261120-geometry-restoration`,
`manifest.json` SHA-256
`03f471d4a4a688ec7d5e90aadb6fae71be006b3d080ab96d52168131556e5a66`, `#292`),
including 61-partition global source acquisition, narrow ring-normalization and
global overlap audit receipts, independent raw-shapefile verification, S3
publication and 1 TB storage headroom verification, remote range probes, and
uncached Pourpoint CAMELS smoke and negative-control checks.

`tdx-catchment-connectivity.json` records the secret-free provenance for the
private TDX Hydro HFX 0.4.0 bundle (`#289`, new immutable prefix
`NGA-TDX-Hydro-20230126-full62-20261121/`) derived on the UZH VM from the
verified `full62-20261027` bundle with single-part `POLYGON` fragment
collection, exhaustively verified zero-area virtual confluence stub
contraction, completeness metadata, S3 readback, and live Pourpoint
end-to-end verification (`cache=False`).
