# TEJ differential splicing assessment

Module authors: Ryan Corbett (@rjcorb)

This module extracts splice events associated with identified TEJs in PBTA cohort, and filters TEJs for those associated with differential splicing events. 

## Usage
### script to run analysis

```
bash run_module.sh
```

## Differential splicing classification strategy

1. Identify rMATS splice events associated with tumor-enriched or oncofetal junctions (TEJs).

   - Events are evaluated only when the TEJ and splice event occur in the same PBTA sample.
   - The workflow considers skipped-exon, retained-intron, alternative 3′ splice-site, and alternative 5′ splice-site events.

2. Compare the tumor event PSI with normal-brain reference PSI matrices from GTEx, postnatal Evo-Devo brain, and normal pediatric brain.

   - Calculate a dPSI for each available normal reference group as `tumor PSI − control PSI`.
   - Use the mean dPSI across available control groups to determine whether the event supports the junction-associated splicing direction.
   - Call inclusion, intron retention, or positive alternative-splice-site events when mean dPSI is greater than 0.1; call exon skipping or negative alternative-splice-site events when mean dPSI is less than −0.1.
   - For events with near-constitutive control PSI, retain directionally supported events based on tumor PSI alone: tumor PSI greater than 0.1 when mean control PSI is below 0.025 for inclusion-like events, or tumor PSI below 0.9 when mean control PSI is above 0.975 for skipping-like events.

3. Resolve multiple differential-splicing events associated with the same sample–junction pair.

   - Retain the sole qualifying event when only one is present.
   - When multiple qualifying events have the same event type, retain the event with the greatest absolute mean dPSI.
   - When qualifying events have different types, retain the event with the greatest absolute mean dPSI.
   - For junctions with no matched control PSI, retain PSI-supported events using the same direction-specific tumor PSI cutoffs and select the event with the greatest directional PSI deviation.

4. Assign a junction-level consensus event type.

   - Across resolved sample-level calls, select the most frequent event type for each junction.
   - Break ties using the largest absolute median dPSI.
   - Only TEJs with a resolved differential-splicing event are retained in the final differential-splicing output.

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
