"""
Parameter register.

Every constant in the generator that can move a reported conclusion lives here
with four things: its value, its unit, where the value came from, and the range
it could plausibly take. Nothing may be load-bearing and unsourced at the same
time.

Provenance codes
  LIT   a published value or range, cited
  PAPER taken from the Temple preprint's own description of its protocol/apparatus
  CAL   set by calibrating the SIMULATED REFERENCE CHANNEL against a statistic the
        preprint reports about its OWN reference instrument (never against the
        device under test)
  ASSUM chosen by judgement, not sourced. Any conclusion that depends on an ASSUM
        parameter must be reported with its robustness range, never as a number.

`lo`/`hi` are the range the robustness harness sweeps. For LIT entries they are
the reported spread; for ASSUM entries they are deliberately wide.
"""

P = {
  # ---- cerebrovascular control -------------------------------------------
  "co2_react_mean": dict(v=0.034, unit="fractional CBF change per mmHg PaCO2",
      lo=0.025, hi=0.045, prov="LIT",
      note="hypercapnic cerebrovascular CO2 reactivity, healthy adults"),
  "co2_react_hypo_ratio": dict(v=0.55, unit="ratio", lo=0.40, hi=0.75, prov="LIT",
      note="hypocapnic reactivity is markedly lower than hypercapnic; asymmetry well described"),
  "autoreg_gain": dict(v=0.0015, unit="fractional MCAv per mmHg CPP", lo=0.0005, hi=0.004,
      prov="ASSUM", note="static autoregulation is effective but not perfect; gain not pinned"),
  "mcav_supine_extra": dict(v=0.060, unit="fraction", lo=0.02, hi=0.10, prov="ASSUM",
      note="supine MCAv rise not explained by CO2 or the modelled CPP term"),
  "mcav_co_coupling": dict(v=0.035, unit="fractional MCAv per unit normalised CO",
      lo=0.01, hi=0.08, prov="ASSUM", note="cerebral flow follows cardiac output weakly"),
  "hydrostatic": dict(v=0.77, unit="mmHg per cm of column", lo=0.74, hi=0.80, prov="LIT",
      note="blood column pressure gradient, standard value"),

  # ---- systemic haemodynamics --------------------------------------------
  # THE HINGE. The postural floor swings from 0.33 to 0.86 across this range.
  "supine_sv_rise": dict(v=0.52, unit="fractional stroke-volume rise, supine vs standing",
      lo=0.40, hi=0.70, prov="LIT",
      note="THE HINGE: its sign decides the postural conclusion. Zaidi 2000 (Heart 83:181, "
           "doi 10.1136/heart.83.2.181, n=20, impedance cardiography): head-up tilt lowers "
           "cardiac output 17-20% and stroke volume 28-34%. Wieling 1998 (Clin Sci 94:347, "
           "doi 10.1042/cs0940347): CO -26+/-10%, SV -39+/-9% after 1 min upright. "
           "Mechanism: Sjostrand 1953, ~640 mL (~11% of blood volume) pools in the legs. "
           "Every source agrees supine CO EXCEEDS standing CO; only the magnitude varies. "
           "v=0.52 reproduces a supine/standing CO ratio of about 1.23."),
  "orthostatic_hr": dict(v=14.0, unit="bpm, standing above supine", lo=8.0, hi=22.0,
      prov="LIT", note="standard orthostatic heart-rate response in healthy young adults"),
  "map_exercise_gain": dict(v=25.0, unit="mmHg at 100% VO2max", lo=15.0, hi=35.0,
      prov="LIT", note="MAP rise across the aerobic range"),
  "co_exercise_gain": dict(v=0.35, unit="fractional stroke-volume rise with workload",
      lo=0.20, hi=0.55, prov="LIT", note="stroke volume plateaus around 40-50% VO2max"),

  # ---- ventilatory control ------------------------------------------------
  "co2_rise_slope": dict(v=7.0, unit="mmHg per unit fractional workload", lo=4.0, hi=10.0,
      prov="ASSUM", note="metabolic CO2 load below the ventilatory threshold"),
  "hypervent_threshold": dict(v=0.52, unit="fraction of VO2max", lo=0.45, hi=0.70,
      prov="LIT", note="ventilatory threshold in untrained-to-active adults"),
  "hypervent_slope": dict(v=55.0, unit="mmHg per unit workload above threshold",
      lo=25.0, hi=80.0, prov="ASSUM", note="sets how far MCAv falls at high intensity"),

  # ---- the optical measurement situation ---------------------------------
  "beta_extracranial": dict(v=0.50, unit="fraction of the temple signal that is extracranial",
      lo=0.25, hi=0.75, prov="ASSUM",
      note="the preprint states the source-detector separation and detector geometry "
           "were not disclosed to the investigators, so this cannot be sourced at all"),
  "ptt_offset": dict(v=0.026, unit="s, extracranial pulse arrival after intracranial",
      lo=0.010, hi=0.045, prov="ASSUM",
      note="load-bearing: with zero offset no filter can separate the two components"),
  "ptt_vasomotor": dict(v=0.022, unit="s of additional delay at full vasoconstriction",
      lo=0.005, hi=0.040, prov="ASSUM", note="pulse transit time moves with vascular tone"),
  "skin_supine_rise": dict(v=0.10, unit="fraction", lo=0.0, hi=0.25, prov="ASSUM",
      note="scalp perfusion rises supine: hydrostatic, plus loss of orthostatic "
           "cutaneous vasoconstriction"),
  "skin_co_coupling": dict(v=0.45, unit="fractional skin flow per unit normalised CO",
      lo=0.20, hi=0.80, prov="LIT",
      note="cutaneous beds are not autoregulated, so they follow perfusion pressure "
           "and cardiac output far more directly than the brain does"),
  "skin_thermal_gain": dict(v=2.4, unit="fractional rise at full thermal load",
      lo=1.2, hi=4.0, prov="LIT", note="cutaneous vasodilation with rising core temperature"),

  # ---- reference instrument ----------------------------------------------
  "tcd_noise": dict(v=0.9, unit="cm/s, additive", lo=0.4, hi=2.0, prov="CAL",
      note="set so simulated MCAv transition magnitudes match those the preprint "
           "reports for its own TCD channel"),
  "tcd_drift_sd": dict(v=0.0015, unit="fractional random-walk step", lo=0.0005, hi=0.004,
      prov="CAL", note="probe locked in a headframe, insonation depth and angle fixed, "
                       "not reacquired within a session"),
  "tcd_dropout": dict(v=0.02, unit="fraction of seconds", lo=0.0, hi=0.06, prov="ASSUM"),

  # ---- protocol (fixed by the preprint, not free) ------------------------
  "tp_cycle_blocks": dict(v=(240, 240, 600), unit="s", lo=None, hi=None, prov="PAPER"),
  "tp_posture_blocks": dict(v=(240, 240, 240, 240), unit="s", lo=None, hi=None, prov="PAPER"),
  "n_analysed": dict(v=20, unit="subjects per protocol", lo=None, hi=None, prov="PAPER"),
  "w_max_mean": dict(v=0.50, unit="fraction of VO2max during the cycling block",
      lo=0.34, hi=0.68, prov="ASSUM",
      note="'maximum self-selected cadence' with no stated resistance; intensity inferred"),
}


def free_parameters():
    """Parameters the robustness harness may sweep."""
    return {k: p for k, p in P.items() if p["lo"] is not None and p["hi"] is not None}


def unsourced():
    return {k: p for k, p in P.items() if p["prov"] == "ASSUM"}


if __name__ == "__main__":
    from collections import Counter
    c = Counter(p["prov"] for p in P.values())
    print(f"{len(P)} registered parameters: " +
          ", ".join(f"{k} {v}" for k, v in sorted(c.items())))
    print(f"\n{len(unsourced())} are unsourced (ASSUM) and may not carry a reported number "
          f"without a robustness range:\n")
    for k, p in unsourced().items():
        print(f"  {k:24s} {str(p['v']):>8}  range {p['lo']}-{p['hi']}")
