# TESJ differential splicing assessment

Module authors: Ryan Corbett (@rjcorb)

This module extracts splice events associated with identified TESJs in PBTA cohort, and filters TESJs for those associated with differential splicing events. 

## Usage
### script to run analysis

```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to run analysis
* `01-get-junction-splice-events.R` extract TESJ-associated splice event PSIs from PBTA rMATS 
* `02-get-junction-diff-splicing.R` calculate TESJ-associated splice event dPSIs, and filter for TESJs associated with differential splicing

## Directory structure
```
.
├── 01-get-junction-splice-events.R
├── 02-get-junction-diff-splicing.R
├── README.md
├── results
│   ├── tumor-enriched-oncofetal-diff-splice-junctions.qs2
│   ├── tumor-enriched-oncofetal-junction-diff-splice-event-annotation.tsv.gz
│   └── tumor-enriched-oncofetal-junction-splice-events.qs2
├── run_module.sh
└── util
    └── rmats-processing-functions.R
```