import json, os, sys, time
import pandas as pd
from family_hits import family_hits


def _done(path):
    if not os.path.exists(path):
        return set()
    with open(path) as fh:
        return {(d["A"], d["B"]) for d in map(json.loads, fh)}


def run_batch(pairs_tsv, out_dir, min_escore=0.4, pause=1.0):
    os.makedirs(out_dir, exist_ok=True)
    S, H, E = (os.path.join(out_dir, f) for f in ("summary.jsonl", "hits.jsonl", "errors.jsonl"))
    pairs = pd.read_csv(pairs_tsv, sep="\t", header=None, usecols=[0, 1], names=["A", "B"])
    pairs = pairs.apply(lambda c: c.str.split("-").str[0]).drop_duplicates()
    done = _done(S)
    todo = [(a, b) for a, b in pairs.itertuples(index=False) if (a, b) not in done]
    print(f"{len(pairs)} pairs, {len(done)} done, {len(todo)} to run")

    with open(S, "a") as fs, open(H, "a") as fh, open(E, "a") as fe:
        for i, (a, b) in enumerate(todo, 1):
            t0 = time.time()
            try:
                summ, hits = family_hits(a, b, min_escore=min_escore)
            except Exception as ex:  # log and move on; retried next run
                fe.write(json.dumps({"A": a, "B": b, "error": repr(ex)}) + "\n"); fe.flush()
                continue
            for rec in hits.to_dict("records"):
                fh.write(json.dumps({"pair_A": a, "pair_B": b, **rec}) + "\n")
            fs.write(json.dumps(summ) + "\n")      # summary last = pair is complete
            fh.flush(); fs.flush()
            if i % 10 == 0 or i == len(todo):
                print(f"  {i}/{len(todo)}  {a}-{b}  {summ['status']}  hits={summ.get('n_hits', 0)}", flush=True)
            time.sleep(max(0, pause - (time.time() - t0)))   # be polite to STRING/InterPro


def load(out_dir):
    """(summary, hits) DataFrames; keeps the last record if a pair was rerun."""
    s = pd.read_json(os.path.join(out_dir, "summary.jsonl"), lines=True).drop_duplicates(["A", "B"], keep="last")
    hp = os.path.join(out_dir, "hits.jsonl")
    h = pd.read_json(hp, lines=True) if os.path.getsize(hp) else pd.DataFrame()
    if len(h):  # drop hits from interrupted pairs that never got a summary line
        h = h.merge(s[["A", "B"]], left_on=["pair_A", "pair_B"], right_on=["A", "B"]).drop(columns=["A", "B"])
        h = h.drop_duplicates()
    return s, h
