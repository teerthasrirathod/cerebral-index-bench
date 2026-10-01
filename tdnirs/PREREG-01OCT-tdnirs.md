# Pre-registration: the bench against Temple's TD-NIRS preprint

Written 1 October 2026, 08:40 IST, before any simulation in this file was run.
Target: Gulati, Rudaeva, Dutta, Rogers, Prajapat, Kumar, Gupta, Goyal, Boas,
"Evaluating Brain Flow Index from a temple-worn wearable against depth-resolved
time-domain NIRS ...", bioRxiv 10.64898/2026.09.17.751664, posted 24 Sep 2026.

Their numbers are known at the time of writing. That is unavoidable. What is fixed
here, before running, is which of their numbers are used to set parameters and
which are held out as predictions. Held-out numbers are never used to tune.

## What is being tested

The bench was built on 17 Sep around one claim about architecture: a temple index
is some function of a signal that mixes an intracranial and an extracranial term.
The question now is narrower and sharper. Given their published protocol and their
published analysis pipeline, what kind of index reproduces their table?

Two candidate architectures, both one-parameter-family members, neither a claim
about their proprietary formula:

- FLOW: the index follows a fast, autoregulation-speed haemodynamic response
  (first-order, tau_f = 6 s), i.e. what a pulsatile, MCAv-like flow signal does.
- VOLUME: the index follows the same slow response the brain-layer HbO shows.

Both pass through a causal trailing window of N seconds (any real-time 1 Hz
device output is causal) plus device noise of SD sigma_d.

## Fixed inputs (provenance)

PAPER  protocol timing: HDT 5 x 180 s, squat 5 x 90 s, supine 5 x 180 s,
       task in blocks 2 and 4. 20 sessions per protocol.
PAPER  analysis: 1 Hz grid; forward-backward EMA alpha 0.35 on both traces;
       Pearson at every integer lag -30..+30 s, keep the maximum; zero-lag Pearson;
       transition windows +/-90 s (HDT, supine), +/-45 s (squat) on z-scored traces;
       Bland-Altman difference = index minus NIRS.
CAL    brain-layer response kinetics from their own description of their own
       reference (Sec 3.1): half of peak within 10-20 s, peak at 60-90 s.
       -> tau_B drawn U(15, 29) s per session. This calibrates the REFERENCE,
       never the device under test.
ASSUM  spontaneous (non-block) brain fluctuation: band-limited noise, tau 10 s.
       Its share of brain variance, 1-f, is unknown, so f is swept over
       {0.60, 0.75, 0.90} and every result is reported per f.
ASSUM  transition ramp 10 s (HDT, supine), 4 s (squat).
ASSUM  NIRS measurement noise: white, SD 0.15 x block amplitude.

## Parameters set from their numbers (calibration targets)

Only two, both from head-down tilt:
  N        set so the median optimal lag in HDT is +14 s (their 13-15 s).
  sigma_d  set so the median lag-adjusted r in HDT is 0.912.

## Held-out predictions (never used to tune)

P1  Squat and supine median lag-adjusted r, against their 0.843 and 0.839.
P2  Zero-lag r in all three, against 0.826 / 0.701 / 0.796.
P3  Ratio zero-lag / lag-adjusted is lowest in squat (shortest blocks).
    Their ratios: 0.906 / 0.832 / 0.949.
P4  Median optimal lag in squat and supine within 13-15 s +/- 3 s.
P5  Bland-Altman standardised bias negative in all three protocols.
    Theirs: -0.511 / -0.441 / -0.251.
P6  Directional concordance >= 95% in all three. Theirs 100 / 100 / 95.

H-ARCH  The architecture test. To produce a +14 s lag, FLOW needs a much longer
        window N than VOLUME, because a fast signal starts out ahead of the slow
        brain response and has to be dragged behind it. Prediction: FLOW with that
        N cannot hold squat (90 s blocks) at or above 0.80 lag-adjusted r, while
        VOLUME can. If FLOW also holds squat >= 0.80, the lag does not
        discriminate the architectures and that is reported as a null.

## The floor, reported as a bracket, no pass/fail

F1  TIMING-ONLY FLOOR. A signal carrying nothing but the protocol timing (the
    task boxcar through a first-order response, best tau in 5-40 s), with the
    fitted device noise, run through their exact pipeline against the brain
    reference. This is what a signal with zero brain content and perfect block
    locking scores.
F2  The same floor with the lag search switched off, to measure how much the
    +/-30 s maximum adds to a no-information signal.
F3  ORACLE CEILING. The true brain latent itself, with device noise, no window.

Headline r is then located between F1 and F3, per protocol, per f.

## What would count as the architecture landing

All of P1, P3, P4, P5, P6 pass, P1 within +/-0.06 of their medians, and H-ARCH
separates the architectures. Anything less is reported as partial, with the
failing prediction named.

## What this cannot do

It cannot recover their formula, their tissue weighting, or anything about the
scalp compartment; the scalp side of their table (0.494 / 0.700 / 0.767) is not
modelled here and no claim is made about it. It treats the TD-NIRS brain layer
as the truth it is compared against, which their own Sec 4.3 says it is not.
