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
* `01-get-junction-splice-events.R` extract junction counts from PBTA rMATS and normalize. 

## Directory structure
```
.
├── 01-get-junction-splice-events.R
├── README.md
├── results
│   └── tumor-enriched-oncofetal-junction-splice-events.qs2
├── run_module.sh
└── util
    └── rmats-processing-functions.R
```