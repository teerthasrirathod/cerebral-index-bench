"""
Per-subject-protocol computation of the candidate metric family.

Nothing here reconstructs a proprietary metric. The three flow candidates and the
three complexity candidates span the plausible construction space, so that the
floor and ceiling can be reported over a family rather than for one guess.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import tbgen as tb

CONDITIONS = ("real", "sham_ic", "sham_symp")
FLOW_CANDIDATES = ("fc_pi", "fc_shape", "fc_ssr")
CX_CANDIDATES = ("cx_wave", "cx_stream", "cx_rri")


def run_cell(sub, pname, pcfg, cond):
    """One subject x protocol x condition. Seeds are fixed per cell, so the three
    conditions differ only in the pathway that condition ablates."""
    dur, kind = pcfg["dur"], pcfg["kind"]
    ph_seed = sub["seed"] + (11 if pname == "cycle" else 29)
    phys = tb.physiology(sub, kind, dur, np.random.default_rng(ph_seed),
                         ablate_symp=(cond == "sham_symp"))
    long_ch, short_ch = tb.sensor(sub, phys, kind,
                                  np.random.default_rng(sub["seed"] + 4441),
                                  ablate_ic=(cond == "sham_ic"))
    ref = tb.references(sub, phys, np.random.default_rng(sub["seed"] + 7717))
    n = len(ref["tcd"])

    fc_pi = tb.flow_candidate(long_ch, mode="pi")[:n]
    return dict(
        phys=phys, raw=long_ch, short=short_ch, ref=ref, n=n,
        fc_pi=fc_pi,
        fc_shape=tb.flow_candidate(long_ch, mode="shape")[:n],
        fc_ssr=tb.flow_candidate(long_ch, short_ch, mode="ssr")[:n],
        cx_wave=tb.entropy_stream(long_ch, tb.FS, 60, "perm")[:n],
        cx_stream=tb.entropy_stream(fc_pi, 1.0, 60, "perm")[:n],
        cx_rri=tb.entropy_beats(long_ch, n, n_beats=140, kind="samp"),
    )
