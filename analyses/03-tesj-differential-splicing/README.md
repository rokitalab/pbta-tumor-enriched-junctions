# TEJ differential splicing assessment

Module authors: Ryan Corbett (@rjcorb)

This module extracts splice events associated with identified TEJs in PBTA cohort, and filters TEJs for those associated with differential splicing events. 

## Usage
### script to run analysis

```
bash run_module.sh
```

## Folder content
* `run_module.sh` shell script to run analysis
* `01-get-junction-splice-events.R` extract PBTA splice events associated with TEJs from rMATS files. 
* `02-get-junction-diff-splicing.R` calculate TESJ-associated splice event dPSIs, and filter for TEJs associated with differential splicing

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