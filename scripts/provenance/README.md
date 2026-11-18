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
