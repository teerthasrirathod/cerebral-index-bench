# cerebral-index-bench

An independent, pre-registered check of a temple-worn cerebral index (Temple's Brain Flow, Temple-BF) against its published validation studies, using only public preprints, open physiological datasets and synthetic data.

Page: https://teerthasrirathod.github.io/cerebral-index-bench/

Author: Teertha Sri Rathod Banoth, undergraduate (BTech Biotechnology, graduating 2027), working independently. ORCID 0009-0004-5193-3035. Not affiliated with, funded by or reviewed by Temple Pvt. Ltd. or any institution.

## Layout

    index.html   the page
    tdnirs/      timing bracket against the TD-NIRS preprint (bioRxiv 10.64898/2026.09.17.751664)
                 PREREG written before the run, FINDINGS after, sim.py reproduces everything in ~90 s
    bench/       synthetic temple-sensor cohort: SPEC.md (written before any result), src/, results/
    co2/         cerebral vs peripheral CO2 response from PhysioNet Cerebral Vasoregulation in Diabetes,
                 and the breath-hold control from OpenNeuro ds003192

## Reproduce

    pip install numpy scipy pandas
    python tdnirs/sim.py                     # sections 1-2 of the page
    python bench/src/generate_cohort.py      # synthetic cohort, ~5.5 min, seeded
    python bench/src/verify.py               # 11 physiology checks; must pass before results are read

The co2/ scripts download PhysioNet records with `wfdb` and need network access.

## Licence

Text CC BY 4.0. Code MIT.
