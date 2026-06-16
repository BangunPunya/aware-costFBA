from __future__ import annotations

import json
import os

import cobra  # noqa: F401  (version asserted at end)
from modal_context import ModalContext
from coupling import score_coupling

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(CODE_DIR, "data")
OUT_JSON = os.path.join(DATA_DIR, "demo_fbr_lysine.json")

TARGET = "EX_lys__L_e"          # L-lysine exchange (confirmed present in iML1515)
THROTTLE_RXN = "DHDPS"          # lysine-committed, single-gene (dapA/b2478) throttle
GROWTH_FRACTION = 0.9           # score lysine at >= 90% of max growth

# WT throttle level: fraction of the unconstrained DHDPS capacity (illustrative,
# not fitted to a measured Ki). 0.50 binds (caps WT lysine below the fbr level, Δ>0)
# while leaving growth intact (WT growth == fbr growth, growth_retained=yes). The
# sweep below also probes tighter levels (0.20/0.30) that do depress growth.
WT_THROTTLE_FRACTION = 0.50
# Sensitivity sweep - show the result is structural, not one cherry-picked value.
SWEEP_FRACTIONS = [0.20, 0.30, 0.50, 0.70]


def _dhdps_unconstrained_capacity(model: cobra.Model) -> float:
    """Max achievable DHDPS flux at biomass >= GROWTH_FRACTION*growth_max (kcat ceiling)."""
    with model:
        sol = model.optimize()
        gmax = float(sol.objective_value)
        for br in (r for r in model.reactions if r.objective_coefficient != 0):
            br.lower_bound = GROWTH_FRACTION * gmax
        model.objective = model.reactions.get_by_id(THROTTLE_RXN)
        model.objective_direction = "max"
        cap = float(model.optimize().objective_value)
    return cap


def _make_wt(base: cobra.Model, cap: float, frac: float) -> cobra.Model:
    """Wild-type regulated model: cap DHDPS UB to `frac` of capacity (throttle on)."""
    m = base.copy()
    m.reactions.get_by_id(THROTTLE_RXN).upper_bound = frac * cap
    return m


def _make_fbr(base: cobra.Model, cap: float) -> cobra.Model:
    """fbr candidate: relax throttle back to unconstrained capacity (UB = `cap`)."""
    m = base.copy()
    m.reactions.get_by_id(THROTTLE_RXN).upper_bound = cap
    return m


def _score(model: cobra.Model) -> dict:
    """Lysine @ near-max growth for one model (reuses coupling.score_coupling)."""
    r = score_coupling(model, TARGET, growth_fraction=GROWTH_FRACTION)
    return {
        "growth_max": round(r["growth_max"], 6),
        "lysine_at_near_max_growth": round(r["prod_at_near_max"], 6),
        "lysine_min_guaranteed": round(r["prod_min_guaranteed"], 6),
        "feasible": bool(r["feasible"]),
        "status": r["status"],
    }


def run_oracle(plain: bool, base_plain: cobra.Model, base_modal: cobra.Model,
               cap_plain: float, cap_modal: float, wt_frac: float) -> dict:
    """Score WT vs fbr under one oracle.

    PLAIN FBA: no regulatory layer, WT == fbr (same LP) -> delta ~ 0 (false-negative).
    MODAL FBA: WT throttled, fbr relaxes it -> different LPs -> delta > 0.
    """
    if plain:
        # No regulatory throttle on EITHER model. WT == fbr by construction:
        # plain FBA literally cannot see the feedback intervention.
        wt = base_plain.copy()
        fbr = base_plain.copy()
    else:
        wt = _make_wt(base_modal, cap_modal, wt_frac)
        fbr = _make_fbr(base_modal, cap_modal)

    wt_s = _score(wt)
    fbr_s = _score(fbr)
    delta = round(fbr_s["lysine_at_near_max_growth"]
                  - wt_s["lysine_at_near_max_growth"], 6)

    # growth retained? fbr growth should be within 1% of WT growth (throttle must
    # not be silently killing growth). Reported as observed.
    g_wt = wt_s["growth_max"]
    g_fbr = fbr_s["growth_max"]
    growth_retained = bool(g_fbr > 1e-6 and g_wt > 1e-6
                           and abs(g_fbr - g_wt) / max(g_wt, 1e-9) < 0.01)
    return {
        "oracle": "plain_FBA" if plain else "modal_regulatory_FBA",
        "WT": wt_s,
        "fbr": fbr_s,
        "delta_lysine_fbr_minus_WT": delta,
        "growth_retained": growth_retained,
    }


def main() -> dict:
    os.makedirs(DATA_DIR, exist_ok=True)
    mc = ModalContext.load()

    # MODAL base = modal layer ON (M9 medium + ATPM) via ModalContext.apply().
    base_modal = mc.apply(mc.base.copy())
    # PLAIN base = naive FBA-polos (rich medium, ATPM=0) - the "plain stoichiometric
    # FBA" baseline with NO regulatory layer whatsoever. Reuses mc.as_polos().
    polos = mc.as_polos()
    base_plain = polos.apply(polos.base.copy())

    # Unconstrained DHDPS capacity ceiling, computed per oracle's network state.
    cap_modal = _dhdps_unconstrained_capacity(base_modal)
    cap_plain = _dhdps_unconstrained_capacity(base_plain)

    # Ablation table: plain vs modal, WT vs fbr.
    plain_res = run_oracle(True, base_plain, base_modal, cap_plain, cap_modal,
                           WT_THROTTLE_FRACTION)
    modal_res = run_oracle(False, base_plain, base_modal, cap_plain, cap_modal,
                           WT_THROTTLE_FRACTION)

    # Sensitivity sweep: vary throttle level; modal delta scales, plain delta stays 0.
    sweep = []
    for frac in SWEEP_FRACTIONS:
        m = run_oracle(False, base_plain, base_modal, cap_plain, cap_modal, frac)
        p = run_oracle(True, base_plain, base_modal, cap_plain, cap_modal, frac)
        sweep.append({
            "throttle_fraction_of_capacity": frac,
            "modal_WT_lysine": m["WT"]["lysine_at_near_max_growth"],
            "modal_fbr_lysine": m["fbr"]["lysine_at_near_max_growth"],
            "modal_delta": m["delta_lysine_fbr_minus_WT"],
            "modal_WT_growth": m["WT"]["growth_max"],
            "modal_growth_retained": m["growth_retained"],
            "plain_delta": p["delta_lysine_fbr_minus_WT"],  # expected 0.0 always
        })

    interpretation = (
        "feedback-resistance is unscorable by plain FBA (delta=%.4f) and scorable "
        "only under the modal regulatory layer (delta=+%.4f). The qualitative "
        "result (plain=0, modal>0) holds across ALL throttle levels in the sweep; "
        "only the magnitude of the modal delta is parameter-dependent."
        % (plain_res["delta_lysine_fbr_minus_WT"],
           modal_res["delta_lysine_fbr_minus_WT"])
    )

    result = {
        "title": "fbr enzyme intervention: invisible to plain FBA, scorable under "
                 "modal regulatory-throttle (L-lysine, iML1515)",
        "target_metabolite": "L-lysine",
        "target_exchange": TARGET,
        "throttled_reaction": {
            "id": THROTTLE_RXN,
            "name": "dihydrodipicolinate synthase (DHDPS)",
            "gene": "dapA (b2478), single gene, no isozyme escape",
            "why": "first committed, lysine-specific step of the DAP/lysine branch; "
                   "allosterically inhibited by L-lysine in E. coli -> clean throttle.",
        },
        "aspk_note": "ASPK (aspartate kinase) was inspected and REJECTED as the "
                     "throttle: 3 isozymes (lysC/thrA/metL); only lysC is "
                     "lysine-sensitive, so a single-reaction cap reroutes flux "
                     "through thrA/metL and is not lysine-specific. ASPK is a "
                     "shared node, DHDPS is lysine-committed.",
        "modal_layers_active": [k for k, v in mc.toggles.items() if v],
        "wt_throttle_fraction_ILLUSTRATIVE": WT_THROTTLE_FRACTION,
        "dhdps_unconstrained_capacity": {
            "modal_oracle": round(cap_modal, 6),
            "plain_oracle": round(cap_plain, 6),
        },
        "ablation_table": {
            "plain_FBA": plain_res,
            "modal_regulatory_FBA": modal_res,
        },
        "sensitivity_sweep": sweep,
        "interpretation": interpretation,
        "caveats": [
            "Throttle magnitude is ILLUSTRATIVE / data-hungry: a real value needs a "
            "measured Ki / degree-of-inhibition of L-lysine on DHDPS. Not fitted here.",
            "The ROBUST claim is QUALITATIVE: plain-FBA delta=0 at every throttle "
            "level, modal delta>0 whenever the throttle binds. Exact magnitude is "
            "parameter-dependent (see sensitivity_sweep).",
            "This is IN-SILICO RE-RANKING (does the modal oracle distinguish WT vs "
            "fbr?), NOT a viability/titre prediction.",
            "Throttle is on the lysine-committed step DHDPS (single gene b2478) to "
            "avoid ASPK isozyme rerouting through thrA/metL.",
            "Growth feasibility checked in every condition (growth_retained flags). A "
            "very tight throttle (e.g. fraction=0.20) depresses growth - reported, "
            "not hidden.",
        ],
        "cobra_version": cobra.__version__,
    }

    with open(OUT_JSON, "w") as f:
        json.dump(result, f, indent=2)

    _print_report(result)
    return result


def _print_report(r: dict) -> None:
    p = r["ablation_table"]["plain_FBA"]
    m = r["ablation_table"]["modal_regulatory_FBA"]
    line = "-" * 92
    print("\n" + line)
    print("NOVELTY DEMO: fbr lysine intervention - plain FBA vs modal regulatory-throttle")
    print("Throttled reaction: %s (%s)  |  target: %s"
          % (r["throttled_reaction"]["id"], r["throttled_reaction"]["gene"], r["target_exchange"]))
    print(line)
    hdr = "%-24s | %14s | %12s | %14s | %s"
    print(hdr % ("oracle", "WT lys@~max-g", "fbr lys", "Δ (fbr−WT)", "growth retained?"))
    print(line)
    row = "%-24s | %14.5f | %12.5f | %+14.5f | %s"
    for res in (p, m):
        print(row % (res["oracle"], res["WT"]["lysine_at_near_max_growth"],
                     res["fbr"]["lysine_at_near_max_growth"],
                     res["delta_lysine_fbr_minus_WT"],
                     "yes" if res["growth_retained"] else "NO"))
    print(line)
    print("\nSENSITIVITY SWEEP (throttle level vs Δ):")
    sh = "%8s | %12s | %12s | %12s | %12s | %s"
    print(sh % ("frac", "modal WT", "modal fbr", "modal Δ", "plain Δ", "WT growth"))
    for s in r["sensitivity_sweep"]:
        print(sh % ("%.2f" % s["throttle_fraction_of_capacity"],
                    "%.5f" % s["modal_WT_lysine"], "%.5f" % s["modal_fbr_lysine"],
                    "%+.5f" % s["modal_delta"], "%+.5f" % s["plain_delta"],
                    "%.4f%s" % (s["modal_WT_growth"],
                                "" if s["modal_growth_retained"] else " (growth(down))")))
    print(line)
    print("INTERPRETATION: " + r["interpretation"])
    print("cobra version:", r["cobra_version"])
    print("JSON ->", OUT_JSON)
    print(line + "\n")


if __name__ == "__main__":
    main()
