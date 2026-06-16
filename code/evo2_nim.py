from __future__ import annotations

import os
import sys
import time

import requests

INVOKE_URL = "https://health.api.nvidia.com/v1/biology/arc/evo2-40b/generate"
STATUS_URL = "https://health.api.nvidia.com/v1/status/"   # + reqid (async poll)
VALID_DNA = set("ACGTN")


KEY_FILE = os.path.expanduser("~/.nvidia_evo2_key")  # fallback if $NVIDIA_API_KEY unset


def _read_key() -> str:
    key = os.environ.get("NVIDIA_API_KEY", "").strip()
    if key:
        return key
    if os.path.exists(KEY_FILE):
        k = open(KEY_FILE).read().strip()
        if k:
            return k
    return ""


def _headers():
    key = _read_key()
    if not key:
        raise RuntimeError(
            "No API key found. Set it via either:\n"
            "    echo 'nvapi-xxxx' > ~/.nvidia_evo2_key\n"
            "    export NVIDIA_API_KEY=nvapi-xxxx\n"
            "Free key: build.nvidia.com/arc/evo2-40b ('Get API Key').")
    return {"Authorization": f"Bearer {key}", "Accept": "application/json",
            "Content-Type": "application/json"}


def generate_dna(sequence: str, num_tokens: int = 100, temperature: float = 0.7,
                 top_k: int = 3, top_p: float = 0.0, random_seed: int | None = None,
                 timeout: int = 300, poll_every: float = 5.0) -> dict:
    """Generate a DNA continuation of `sequence`. Returns the Evo2 response dict."""
    seq = sequence.strip().upper()
    bad = set(seq) - VALID_DNA
    if bad:
        raise ValueError(f"prompt contains non-DNA characters: {bad} (only A/C/G/T/N)")
    payload = {"sequence": seq, "num_tokens": int(num_tokens),
               "temperature": float(temperature), "top_k": int(top_k),
               "top_p": float(top_p), "enable_sampled_probs": False,
               "enable_logits": False}
    if random_seed is not None:
        payload["random_seed"] = int(random_seed)

    sess = requests.Session()
    r = sess.post(INVOKE_URL, headers=_headers(), json=payload, timeout=120)

    # async: 202 + poll status endpoint
    deadline = time.time() + timeout
    while r.status_code == 202:
        reqid = r.headers.get("NVCF-REQID") or r.headers.get("nvcf-reqid")
        if not reqid:
            break
        if time.time() > deadline:
            raise TimeoutError("Evo2 NIM poll timeout")
        time.sleep(poll_every)
        r = sess.get(STATUS_URL + reqid, headers=_headers(), timeout=120)

    if r.status_code != 200:
        raise RuntimeError(f"Evo2 NIM HTTP {r.status_code}: {r.text[:300]}")
    return r.json()


def _self_check():
    """Validate setup without calling the API (check key)."""
    try:
        _headers()
        print("[ok] NVIDIA_API_KEY read")
    except RuntimeError as e:
        print(f"[!] {e}")
        return False
    return True


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        print("\nCheck setup:  python3 evo2_nim.py --check")
        sys.exit(0)
    if sys.argv[1] == "--check":
        sys.exit(0 if _self_check() else 1)

    prompt = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    print(f"[evo2-40b] prompt={prompt[:40]}... num_tokens={n}", flush=True)
    out = generate_dna(prompt, num_tokens=n)
    gen = out.get("sequence", "")
    cont = gen[len(prompt):] if gen.startswith(prompt) else gen
    print(f"[evo2-40b] elapsed_ms={out.get('elapsed_ms')}")
    print(f"[evo2-40b] generated ({len(cont)} nt continuation):")
    print(cont)
    # write FASTA
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "evo2_output.fasta"), "w") as f:
        f.write(f">evo2_40b_gen num_tokens={n}\n{gen}\n")
    print("[evo2-40b] written -> evo2_output.fasta")
