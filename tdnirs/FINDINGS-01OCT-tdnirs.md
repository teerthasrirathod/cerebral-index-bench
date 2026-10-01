# The bench against Temple's TD-NIRS preprint: what landed

1 October 2026. Rules in `tdnirs/PREREG-01OCT-tdnirs.md`, written before the run.
Code `tdnirs/sim.py`, raw output `tdnirs/results_tdnirs.json`. Runtime 86 s.
Target: bioRxiv 10.64898/2026.09.17.751664 (Gulati ... Boas, posted 24 Sep 2026).

Calibrated on head-down tilt only (two numbers: lag +14 s, r 0.912). Everything
below for squat and supine, and every secondary statistic, is a held-out prediction.
f = share of brain-layer variance that is block-locked; unknown, so swept.

---

## 1. The architecture finding: a ~30 s trailing window

To put Temple-BF 14 s behind the brain layer, a VOLUME-type index needs a causal
trailing window of **30 s, at every f**. A FLOW-type index needs 58-70 s.

One parameter, fitted on one protocol's lag. What it predicts elsewhere:

| prediction | bench (VOLUME, f 0.75) | theirs | verdict |
|---|---|---|---|
| P4 lag, squat / supine | +14.5 / +14.0 s | 13-15 s | **pass** |
| P3 zero-lag/lagged ratio lowest in squat | 0.911 / **0.872** / 0.915 | 0.906 / **0.832** / 0.949 | **pass** |
| P2 zero-lag r, HDT | 0.830 | 0.826 | **pass** |
| P5 Bland-Altman bias negative, all three | -0.37 / -0.66 / -0.34 | -0.51 / -0.44 / -0.25 | **pass** on sign; magnitude order wrong |
| P6 concordance >= 95% | 100 / 95 / 100 | 100 / 100 / 95 | **pass** |
| P1 lag-adjusted r, squat / supine | 0.906 / 0.909 | 0.843 / 0.839 | **fail**, over by 0.06-0.07 |
| P2 zero-lag r, squat / supine | 0.790 / 0.832 | 0.701 / 0.796 | **fail** on squat |

So the timing structure of their table lands: constant lag across protocols, the
largest zero-lag penalty in the shortest blocks, an index whose standardized
transition response is smaller than the reference's. All of that follows from a
single trailing window of about 30 s. **Read plainly: Temple-BF behaves like a
signal averaged over roughly the last half-minute.** That is an inference about
their output's timing, not their formula.

What does not land: the bench is too good in squat and supine. With the HDT noise
level it predicts ~0.91 everywhere; they get 0.84. Something degrades Temple-BF
in those two protocols that timing alone does not. Exploratory and post hoc, not
pre-registered: raising device noise 1.4-1.6x in squat and supine brings
lag-adjusted r to their values and then predicts zero-lag r of 0.72-0.74 (squat,
theirs 0.701) and 0.78-0.80 (supine, theirs 0.796) without further tuning. That
points at protocol-specific degradation, plausibly motion or a superficial
contribution, but the bench cannot say which. That is a question for their data.

## 2. H-ARCH: the lag does not separate flow from volume

Pre-registered expectation: FLOW, needing a ~60 s window, would collapse on the
90 s squat blocks. It does not: 0.894-0.905. **Null, as written down in advance.**
The lag alone cannot say whether Temple-BF is flow-like or volume-like. Not
pre-registered, and noted only: FLOW overshoots their HDT zero-lag r (0.867-0.874
against 0.826), VOLUME at f 0.75 hits it (0.830). Weak, one number.

## 3. The bracket

Floor: a signal containing only the protocol timing, no brain, same device noise,
their exact pipeline, best response time constant chosen for it. Ceiling: the true
brain latent with the same device noise, no window.

| | floor f 0.60 / 0.75 / 0.90 | ceiling | theirs | position at f 0.75 |
|---|---|---|---|---|
| HDT | 0.750 / 0.817 / 0.879 | 0.92-0.94 | **0.912** | 0.88 of the way up |
| squat | 0.756 / 0.825 / 0.886 | 0.91-0.93 | **0.843** | 0.18 |
| supine | 0.750 / 0.817 / 0.879 | 0.92-0.94 | **0.839** | 0.20 |

- **Head-down tilt is the only protocol where Temple-BF sits clearly above what
  protocol timing alone produces**, at every f.
- In squat and supine it sits just above the timing-only floor at f 0.75, and
  below it at f 0.90. If the device is noisier there (section 1), the floor drops
  with it, so the squat/supine verdict is conditional on that.
- This matches their own brain-minus-scalp result protocol for protocol: large in
  HDT (0.333), small in squat (0.110), not significant in supine (0.036). Two
  independent routes, their scalp contrast and this timing floor, agree on which
  protocol carries the evidence.
- The timing-only floor's best lag is negative (-4 to -10 s): a pure timing
  signal runs ahead of the slow brain response. Temple-BF's +14 s is therefore
  not a property of the protocol; it is the device's averaging.

## 4. Consistency with the earlier work

The synthetic cohort in `bench/` was built on 17 September against the
transcranial Doppler preprint. Its robust conclusions under standing/supine are
sign claims only (the bench's own SPEC forbids quoting its floor magnitudes):
posture moves cerebral and extracranial perfusion together, so a postural
challenge cannot attribute a signal to the brain. The TD-NIRS supine result
(brain vs scalp p = 0.082) and the floor above say the same thing from different
data. Nothing from the earlier work is contradicted.

## 5. What their data would settle, and nothing else can

1. f, the block-locked share of brain-layer variance. Every floor above depends
   on it, and their session traces give it directly.
2. The residual correlation: Temple-BF against brain-layer HbO after the
   block-locked mean is removed. That is the test of brain tracking beyond the
   protocol, and it is the one number a timing-only signal cannot fake.
3. Whether the squat/supine shortfall is device noise (motion) or content.
4. The window: a step response of Temple-BF v10.0 to any input would confirm or
   kill the ~30 s figure in one plot.

## 6. Limits

1 Hz reduced model, not the 64 Hz generator. Brain-layer HbO is treated as truth,
which their Sec 4.3 says it is not. The scalp side of their table is not modelled.
Spontaneous brain fluctuation and NIRS noise are assumed. The two architectures
are family members, not their algorithm.
