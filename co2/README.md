# CO2: cerebral and peripheral responses move in opposite directions

Data: PhysioNet Cerebral Vasoregulation in Diabetes 1.0.0, Day-2 head-up-tilt records, 31 subjects aged 55-75.
Breath-hold control: OpenNeuro ds003192 (EuskalIBUR) physiological recordings, 7 subjects.

| | response to end-tidal CO2 | consistency |
|---|---|---|
| cerebral MCA velocity (MFV) | +2.3 to +3.3 %/mmHg over 360 analysis settings, median 3.02 | positive in 21 of 21 |
| peripheral artery velocity (MFV) | -1.04 %/mmHg | negative in 16 of 21, sign-flip permutation p = 0.011 |
| arterial pressure | R2 0.000 | none |

Heart rate adds 0.0003 of R2 to the cerebral model; arterial pressure adds 0.0001.
Peripheral effect in the arithmetic epoch mean: -0.14 %/mmHg (R2 0.005), because it lives in diastole.
Breath-hold: O2 falls with CO2 in 32 of 41 runs (p = 0.0002), so breath-hold is not a single-variable CO2 stimulus.
Optical arm: two attempts from open data (prefrontal NIRS under orthostatic challenge; scanner pulse oximetry under breath-hold) returned nulls.

Velocities, not flows; relative change only.

## Files

    vasoreg/reduce.py      every reduction, each declared in one place; the time base is applied once, here
    vasoreg/run_all.py     cohort run; vasoreg/sweep.py the 360-setting sensitivity sweep
    vasoreg/analyze_cohort.py, tail_test.py
    vasoreg/epochs_all.jsonl, sweep_raw.jsonl, sweep_fits.json   recorded outputs
    breathhold/bhreduce.py, run_bh.py, bh_rows.jsonl.gz          breath-wise features, no task timing used
