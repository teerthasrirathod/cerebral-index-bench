# Synthetic temple-wearable cohort: generator specification

Written before any result is computed. Rules stated here are binding on the analysis;
anything changed afterwards is recorded as a change, not silently edited.

## 0. Purpose

Temple Pvt. Ltd. reports two numbers that cannot be checked from outside:

- **Brain Flow**: median within-subject correlation ~0.80 against transcranial Doppler
  middle cerebral artery velocity (MCAv), n=23 healthy adults, 20 analysed per protocol,
  two protocols (stationary cycling; standing-to-supine postural shift). Preprint exists.
- **Entropy**: r = 0.93 against metabolic-cart VO2 across >100 cardio sessions, versus
  r = 0.55 for heart rate. Scale 1-250, updated every second. No methods published.

Neither number means anything against zero. Both need a two-sided bracket: a floor
(what the number reaches with the cerebral/metabolic information removed but the
hardware and the confounds intact) and a ceiling (the best any estimator could reach
given the reference standard's own noise and the stream rate). A synthetic cohort with
known ground truth is the only way to build that bracket without their data.

This generator is not a model of Temple's device. It is a model of the measurement
situation their device is in.

## 0.1 What can be said about a metric whose formula is not public

Neither Brain Flow nor Entropy is disclosed, and nothing here reconstructs either.
The bracket does not need them.

Whatever function f a device applies, f can only use information present in the
signal. That gives two algorithm-independent bounds:

- **The floor is algorithm-independent.** Ablate the cerebral component of the
  sensor stream, hold the scalp, heart rate, motion, drift and noise bit-identical,
  and hand the ablated stream to f. Whatever f scores there is what f scores with
  no brain in it. Bounding a reported number therefore does not require knowing f;
  it requires showing that a *spanning family* of plausible f reaches that number
  on an ablated stream. If a whole family does, the burden moves: the metric's
  owner must show theirs does something no member of the family does.
- **The ceiling is algorithm-independent.** It is set by the reference standard's
  own noise, its dropout and the stream rate. No f beats it.

The middle of the bracket is where the algorithm lives, and it is not filled here.
That is the point: the bracket is honest because it brackets *over* algorithms
instead of guessing one.

One further handle exists for Entropy specifically. The public description - a
complexity measure, 1 to 250, updated every second, from a temple sensor - does not
fix the formula but does narrow the *input*. A per-second complexity measure on this
hardware is computed on one of three things: the raw waveform, the beat-interval
series recovered from it, or a derived low-rate stream. Those three are implemented
here as `cx_wave`, `cx_rri` and `cx_stream`. *Entropy of what* is a much smaller
question than *entropy how*, and the three inputs behave so differently that the
question answers part of itself.

Accordingly nothing in the code is named after their metrics. The flow candidates
are `fc_pi`, `fc_shape`, `fc_ssr`; the complexity candidates are `cx_wave`,
`cx_stream`, `cx_rri`. Each is one member of a family, not a claim about anyone's
product.

## 1. Cohort

24 synthetic subjects. Each completes both protocols. Subject-level parameters drawn
once and held fixed across protocols.

| Parameter | Distribution | Anchor |
|---|---|---|
| age | U(22, 45) yr | healthy-adult sample |
| body mass | N(68, 11) kg | |
| HR_rest (supine) | N(58, 7) bpm | |
| HR_max | 208 - 0.7*age, +/- N(0,5) | Tanaka |
| VO2_rest | 3.5 ml/kg/min fixed | 1 MET |
| VO2_max | N(38, 7) ml/kg/min, clipped [24, 58] | untrained-to-active |
| MCAv_rest | N(58, 9) cm/s, clipped [40, 80] | adult MCA resting range |
| CO2 reactivity | N(3.4, 0.6) %/mmHg hypercapnic; 0.55x that hypocapnic | asymmetric reactivity |
| beta_extracranial | Beta-shaped on [0.25, 0.75], mean ~0.45 | scalp/superficial temporal fraction of a temple-surface signal |
| head height above heart, standing | N(28, 3) cm | hydrostatic column |
| probe gain error (TCD) | lognormal, CV 12% | insonation angle, fixed per subject per protocol |

`beta_extracranial` is the parameter the whole exercise turns on. A sensor on the temple
sits over the superficial temporal artery and scalp. Whatever fraction of its signal is
extracranial carries heart rate, skin blood flow and thermoregulation, and carries no
cerebral information at all.

## 2. Protocols

### A. Cycling ramp, 16 min, 64 Hz + 1 Hz

    0-180 s     seated rest
    180-360 s   30% VO2max
    360-600 s   50% VO2max
    600-840 s   70% VO2max
    840-960 s   recovery, unloaded

Cadence 80 rpm (1.333 Hz) with subject-level jitter, present from 180 s.

### B. Postural shift, 10 min, 64 Hz + 1 Hz

    0-300 s     quiet standing
    300-310 s   transition
    310-600 s   supine

## 3. Latent physiology

Workload W(t) in [0,1] of VO2max.

**VO2**: VO2(t) = VO2_rest + (VO2_max - VO2_rest) * (W conv mono-exponential, tau = 35 s,
cardiodynamic delay 15 s). Phase II kinetics; tau drawn N(35, 6) s per subject.
Recovery tau 42 s. Supine VO2 falls ~8% below standing.

**HR**: HR(t) = HR_rest + (HR_max - HR_rest) * (W conv exp, tau = 28 s)
         + cardiac drift 0.04 bpm/s accumulated above 40% VO2max
         + respiratory sinus arrhythmia + N(0, 1.2) bpm beat noise.
Standing adds +14 bpm over supine (orthostatic). At the supine transition HR falls with
tau = 20 s.

**PaCO2**: +2.5 mmHg at moderate intensity (metabolic CO2 load), then falls as
hyperventilation begins above ~60% VO2max, reaching -7 mmHg at 70% VO2max.
This is what makes MCAv non-monotone in workload.

**MAP**: rises 0.25 mmHg per % VO2max. Supine transition: cerebral perfusion pressure
rises by rho*g*h, i.e. 0.77 mmHg per cm of head-above-heart height that is removed.

**MCAv** (ground truth):

    MCAv(t) = MCAv_rest * (1 + k_CO2 * dPaCO2(t)) * (1 + 0.0028 * dCPP(t)) * (1 + 0.15 * dCO_norm(t))

with a 6 s autoregulatory low-pass on the CPP term. Consequences that are *not* fitted
and must show up on their own:

- an inverted-U across exercise intensity: MCAv up ~12-18% at moderate, back toward or
  below baseline at 70% VO2max
- supine MCAv ~10-14% above standing, while HR falls. This is the dissociation.

## 4. Forward sensor model, 64 Hz

    s(t) = (1-beta) * A_ic(t) * P_ic(phi(t))      intracranial-coupled pulsatile
         +   beta   * A_ec(t) * P_ec(phi(t))      extracranial: scalp + superficial temporal
         +          m(t)                          motion
         +          d(t)                          baseline drift
         +          eps(t)                        sensor noise

- phi(t): cardiac phase integrated from instantaneous HR, so both pulsatile terms are
  phase-locked to the heart.
- P_ic: pulse shape with a systolic upstroke and a dicrotic notch whose position moves
  with MAP and vascular tone. A_ic proportional to MCAv(t) with a pulsatility term.
- P_ec: broader, later-peaking scalp waveform. A_ec proportional to skin blood flow,
  which follows core temperature (rises through exercise, keeps rising into recovery)
  and falls on lying down. A_ec carries no cerebral information by construction.
- m(t): during cycling, cadence fundamental plus 2nd and 3rd harmonics, amplitude
  scaling with workload; during the postural shift, one large transient at 300-310 s.
- d(t): random walk plus a slow thermal ramp (sensor contact and skin temperature).
- eps(t): white, SNR set per subject.

## 5. Device outputs

Two channels are exposed, because any device that separates cerebral from scalp
content needs both:

- **long channel**, 64 Hz: the mixed signal above.
- **short-separation channel**, 64 Hz: scalp and superficial temporal artery only,
  sharing the motion and most of the drift, with independent sensor noise.

The extracranial pulse arrives later than the intracranial one and the offset moves
with vasomotor tone (26-48 ms). Without that, the two components sit at identical
cardiac phase and no filter could separate them - which would make short-separation
regression impossible by construction and the whole exercise rigged.

**Flow candidates**, 1 Hz, all normalised to the subject's first 60 s:

- `fc_pi`    perfusion index: beat pulse amplitude over DC light level. The standard
             AC/DC construction, with no separation of cerebral from scalp content.
- `fc_shape` heart-rate-invariant morphology: beats phase-normalised onto a fixed
             grid and ensemble-averaged, then early-cycle area over mid-cycle area.
- `fc_ssr`   short-separation-regressed: a 21-tap FIR fitted from the short channel
             to the long channel over the whole record and subtracted. Fitted once
             rather than per window, so it removes the stationary scalp contribution
             and leaves what changes differently from the scalp.

**Complexity candidates**, 1 Hz, each also mapped onto the reported 1-250 scale:

- `cx_wave`   permutation entropy of the 64 Hz waveform, 60 s window.
- `cx_stream` permutation entropy of the 1 Hz flow stream, 60 s window.
- `cx_rri`    sample entropy of the last 140 detected beat intervals. Beat detection
              runs on the noisy trace, so motion corrupts the intervals exactly when
              workload is highest. That corruption is left in.

## 6. Reference channels

- **TCD MCAv**, 1 Hz beat-averaged: true MCAv * per-subject probe gain, * slow insonation
  drift (multiplicative random walk, SD 3%), + N(0, 1.8) cm/s, with 2% of seconds dropped.
- **Metabolic cart VO2**, breath-by-breath: breaths generated at the ventilatory rate,
  each breath's VO2 = true VO2 * (1 + N(0, 0.09)). Breath noise of this size is why the
  averaging window is a free parameter and not a detail.
- **ECG HR**, 1 Hz: true HR + N(0, 0.6) bpm.

## 7. The bracket

**Floor (the TTX analogue).** Two ablations, because the two claims have different
floors.

- `sham_ic`: the intracranial term is replaced by an extracranial-shaped term of the
  same amplitude budget. Same cardiac phase, same scalp waveform, same motion, same
  drift, same noise seed. Zero cerebral information, identical hardware. This is the
  floor for the flow candidates.
- `sham_symp`: the sympathetic drive that suppresses beat-interval complexity is held
  at a constant resting value. Mean heart rate, cerebral flow, motion and noise are
  unchanged; only the pathway by which a complexity measure could carry information
  heart rate does not is removed. This is the floor for the complexity candidates,
  and it is the right one, because the question there is never "is there brain in it"
  but "is there anything in it beyond heart rate".

Whatever correlation survives an ablation is what the candidate scores with the
signal removed.

**Ceiling (the InactiveChannels analogue).** Give an oracle the true latent MCAv (or true
VO2) and degrade it only by what no estimator can undo: the reference channel's own
noise, the dropout, and the stream rate. Correlate oracle against reference. Nothing
measured through that reference beats this.

A reported number is then located as a position in [floor, ceiling], not as a distance
from zero.

## 8. Analysis rules, fixed in advance

- Correlation is **within-subject, within-protocol**, Pearson on the 1 Hz series, then
  the **median across subjects**. This matches what the preprint reports.
- Pooled-across-subject correlation is computed as well, but only as a named contrast,
  never as the headline.
- First 30 s of every protocol dropped (sensor settling).
- The postural transition window 295-315 s is excluded from protocol B correlations and
  analysed separately.
- A subject-protocol is excluded if >10% of reference seconds are missing, or if the
  synthetic TCD trace contains a gap longer than 5 s.
- Entropy's VO2 model: OLS of VO2 on HR, then HR + Entropy. Report delta R-squared.
  Same epochs, same exclusions, for both.
- Lag: no lag correction in the primary analysis. Lag-optimised correlation is reported
  only in the free-parameter sweep, with the chosen lag stated.
- Nothing is reported that has not been computed on both the 1 Hz and the 64 Hz stream.

## 9. What this cannot do

It cannot say what Temple's device measures, what their Brain Flow algorithm is, or
whether their 0.80 is right. It gives the range a number like 0.80 has to be read
against, under a measurement situation built from published physiology rather than
from their data.

## 10. Sweep method: one-at-a-time overstates robustness

The first robustness pass varied one registered parameter at a time between its
range endpoints: 39 settings over 19 parameters. That probes the axes through one
corner of a 19-dimensional space, and it holds coupled parameters at baseline
while moving their partners, which constructs phenotypes that cannot exist. A
subject with a large orthostatic stroke-volume drop necessarily also has a large
compensatory tachycardia and strong standing cutaneous vasoconstriction; they are
one reflex.

The second pass draws 400 coherent joint phenotypes: parameters are grouped into
clusters (orthostatic reflex magnitude, ventilatory responsiveness, exercise
intensity and cardiovascular gain, recording quality), each cluster takes one
latent, and members map monotonically onto their registered ranges from it.

One-at-a-time overstated robustness on every claim it was applied to:

| claim | one-at-a-time | correlated joint |
|---|---|---|
| exercise scalp floor reaches the reported value | 87% (34/39) | 47% |
| exercise heart rate reaches it | 90% (35/39) | 62% |
| posture scalp floor reaches it | 82% (32/39) | 59% |
| posture cardiac output reaches it | 100% (39/39) | 95% |
| posture cardiac-output confound is positive | 100% | 100% |
| posture heart rate is negative | 100% | 100% |

No floor magnitude may be reported as a number. Only the two sign claims survive,
and they survive with margin: across 400 draws the postural cardiac-output
correlation never falls below +0.664 and never changes sign, and heart rate never
rises above -0.585.
