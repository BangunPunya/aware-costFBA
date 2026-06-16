from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Iterable, Optional

DEFAULT_MODEL = "facebook/esm2_t6_8M_UR50D"

# Default location of the dedicated torch env's python (override with ESM2_PYTHON).
_DEFAULT_ESM_PYTHONS = [
    "/home/mp2mk/miniconda3/envs/esm/bin/python",
]

# Recommended fold-safety threshold (delta logP), in WT-MARGINAL mode. A proposal
# PASSES the gate when its score >= THRESHOLD_T. Calibrated in test_esm2_gate.py:
# this T passes all four documented lysC fbr mutants (worst S345F ~= -5.7) while
# rejecting catastrophic backbone/motif breakers (P-loop Gly->Pro/Trp ~= -6.5..-9.8).
# This is a PERMISSIVE guardrail that rejects fold-destroying outliers, NOT a
# discriminative fbr classifier (see esm2_gate_eval.json "interpretation").
THRESHOLD_T = -6.0
DEFAULT_MODE = "wt"


# Proposal container (light; works with operator_mutate.MutationProposal too).
@dataclass
class GateProposal:
    """Minimal mutation descriptor: wt_aa @ position -> mt_aa."""

    position: int  # 1-based residue position
    wt_aa: str
    mt_aa: str

    def key(self) -> str:
        return f"{self.wt_aa}{self.position}{self.mt_aa}"


def _coerce(proposals: Iterable) -> list[GateProposal]:
    """Accept GateProposal, dicts, (pos,wt,mt) tuples, or operator_mutate MutationProposal objects."""
    out: list[GateProposal] = []
    for p in proposals:
        if isinstance(p, GateProposal):
            out.append(p)
        elif isinstance(p, dict):
            out.append(GateProposal(int(p["position"]), p["wt_aa"], p["mt_aa"]))
        elif isinstance(p, (tuple, list)) and len(p) == 3:
            out.append(GateProposal(int(p[0]), str(p[1]), str(p[2])))
        elif hasattr(p, "candidate_mutations"):  # operator_mutate.MutationProposal
            for mt in p.candidate_mutations:
                out.append(GateProposal(int(p.position), p.wt_aa, mt))
        elif hasattr(p, "position") and hasattr(p, "mt_aa"):
            out.append(GateProposal(int(p.position), p.wt_aa, p.mt_aa))
        else:
            raise TypeError(f"Unrecognized proposal: {p!r}")
    return out


# In-process scoring (only runs when torch is importable here).
_MODEL_CACHE: dict = {}


def _load_model(model_name: str):
    if model_name in _MODEL_CACHE:
        return _MODEL_CACHE[model_name]
    from transformers import AutoModelForMaskedLM, AutoTokenizer  # noqa: PLC0415

    tok = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForMaskedLM.from_pretrained(model_name)
    model.eval()
    _MODEL_CACHE[model_name] = (tok, model)
    return tok, model


def _score_inprocess(
    seq: str, items: list[GateProposal], model_name: str, mode: str = "masked"
) -> dict[str, float]:
    """Compute delta-logP scores in the current interpreter (torch present). ESM tokenization is [CLS] r1..rN [EOS], so residue i (1-based) is at token index i. mode='masked' masks each position (one pass per distinct pos); mode='wt' uses a single forward pass over the intact sequence (standard ESM variant-effect score, less noisy, recommended for fold-gating)."""
    import torch  # noqa: PLC0415

    tok, model = _load_model(model_name)
    enc = tok(seq, return_tensors="pt")
    input_ids = enc["input_ids"]

    out: dict[str, float] = {}
    with torch.no_grad():
        if mode == "wt":
            logits_all = model(input_ids).logits[0]  # [L+2, vocab]
            logp_all = torch.log_softmax(logits_all, dim=-1)
            for it in items:
                idx = it.position
                wt = it.wt_aa or seq[it.position - 1]
                wt_id = tok.convert_tokens_to_ids(wt)
                mt_id = tok.convert_tokens_to_ids(it.mt_aa)
                delta = float(logp_all[idx, mt_id] - logp_all[idx, wt_id])
                out[it.key()] = round(delta, 4)
            return out

        by_pos: dict[int, list[GateProposal]] = {}
        for it in items:
            by_pos.setdefault(it.position, []).append(it)
        for pos, group in by_pos.items():
            idx = pos  # CLS offset == 1, so token index == 1-based residue pos
            wt_in_seq = seq[pos - 1]
            masked = input_ids.clone()
            masked[0, idx] = tok.mask_token_id
            logits = model(masked).logits[0, idx]
            logp = torch.log_softmax(logits, dim=-1)
            for it in group:
                wt = it.wt_aa or wt_in_seq
                wt_id = tok.convert_tokens_to_ids(wt)
                mt_id = tok.convert_tokens_to_ids(it.mt_aa)
                delta = float(logp[mt_id] - logp[wt_id])
                out[it.key()] = round(delta, 4)
    return out


# Env discovery + subprocess worker.
def _torch_here() -> bool:
    try:
        import torch  # noqa: F401, PLC0415
        import transformers  # noqa: F401, PLC0415

        return True
    except Exception:  # noqa: BLE001
        return False


def _find_esm_python() -> Optional[str]:
    """Locate a python interpreter that has torch+transformers (the esm env)."""
    cand = os.environ.get("ESM2_PYTHON")
    cands = ([cand] if cand else []) + _DEFAULT_ESM_PYTHONS
    extra = shutil.which("python3.11") or shutil.which("python3.12")
    if extra:
        cands.append(extra)
    for c in cands:
        if c and os.path.exists(c):
            return c
    return None


def _clean_env() -> dict:
    """Env for the worker subprocess: strip Windows PATH WSL may leak (breaks conda shims), keep HF cache on the WSL side."""
    env = dict(os.environ)
    env["PATH"] = "/usr/bin:/bin:" + env.get("PATH", "")
    env.setdefault("HF_HOME", os.path.expanduser("~/.cache/huggingface"))
    env.setdefault("TOKENIZERS_PARALLELISM", "false")
    return env


def _score_subprocess(
    seq: str, items: list[GateProposal], model_name: str, mode: str = "masked"
) -> Optional[dict[str, float]]:
    py = _find_esm_python()
    if not py:
        return None
    payload = {
        "seq": seq,
        "items": [[it.position, it.wt_aa, it.mt_aa] for it in items],
        "model": model_name,
        "mode": mode,
    }
    try:
        proc = subprocess.run(
            [py, os.path.abspath(__file__), "--worker"],
            input=json.dumps(payload),
            capture_output=True,
            text=True,
            env=_clean_env(),
            timeout=600,
        )
    except Exception:  # noqa: BLE001
        return None
    if proc.returncode != 0:
        return None
    try:
        resp = json.loads(proc.stdout.strip().splitlines()[-1])
    except Exception:  # noqa: BLE001
        return None
    if resp.get("ok"):
        return {k: float(v) for k, v in resp["scores"].items()}
    return None


# Public API.
def esm2_available() -> bool:
    """True if scores can be produced here or via a discovered torch env."""
    return _torch_here() or (_find_esm_python() is not None)


def score_proposals(
    seq: str,
    proposals: Iterable,
    model_name: str = DEFAULT_MODEL,
    mode: str = "masked",
) -> Optional[dict[str, float]]:
    """delta-logP for each proposal. mode='masked' (masked-marginal) or 'wt' (wt-marginal, less noisy). Returns {"<wt><pos><mt>": delta} or None if no torch env reachable (DEFERRED)."""
    items = _coerce(proposals)
    if not items:
        return {}
    if _torch_here():
        return _score_inprocess(seq, items, model_name, mode)
    return _score_subprocess(seq, items, model_name, mode)


def masked_marginal(
    seq: str,
    position: int,
    wt_aa: str,
    mt_aa: str,
    model_name: str = DEFAULT_MODEL,
) -> Optional[float]:
    """logP(mt | masked) - logP(wt | masked) at a single position. None if DEFERRED."""
    res = score_proposals(seq, [GateProposal(position, wt_aa, mt_aa)], model_name)
    if res is None:
        return None
    return next(iter(res.values()))


def fold_gate(
    seq: str,
    proposals: Iterable,
    threshold: float = THRESHOLD_T,
    model_name: str = DEFAULT_MODEL,
    mode: str = DEFAULT_MODE,
) -> Optional[list[dict]]:
    """Annotate each proposal with delta-logP score and pass/fail fold-safety verdict (score >= threshold => pass); None if scoring is DEFERRED."""
    items = _coerce(proposals)
    scores = score_proposals(seq, items, model_name, mode)
    if scores is None:
        return None
    out: list[dict] = []
    for it in items:
        s = scores.get(it.key())
        out.append(
            {
                "mutation": it.key(),
                "position": it.position,
                "wt_aa": it.wt_aa,
                "mt_aa": it.mt_aa,
                "score": s,
                "threshold": threshold,
                "fold_pass": (s is not None and s >= threshold),
            }
        )
    return out


# Worker entrypoint (re-exec'd under the esm python).
def _worker_main() -> int:
    raw = sys.stdin.read()
    try:
        req = json.loads(raw)
        seq = req["seq"]
        items = [GateProposal(int(p), w, m) for (p, w, m) in req["items"]]
        model_name = req.get("model", DEFAULT_MODEL)
        mode = req.get("mode", "masked")
        scores = _score_inprocess(seq, items, model_name, mode)
        print(json.dumps({"ok": True, "scores": scores}))
        return 0
    except Exception as exc:  # noqa: BLE001
        print(json.dumps({"ok": False, "error": repr(exc)}))
        return 1


if __name__ == "__main__":
    if "--worker" in sys.argv:
        raise SystemExit(_worker_main())
    # Tiny self-demo
    s = "MSEIVVSKFGGTSVADFDAMNRSADIVLSDANVRLVVLSA"
    print("esm2_available:", esm2_available())
    g = fold_gate(s, [GateProposal(20, "M", "I"), GateProposal(8, "K", "P")])
    print(json.dumps(g, indent=2))
