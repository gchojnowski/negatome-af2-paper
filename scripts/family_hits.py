"""Interacting members of the Pfam families of two UniProt accessions (STRING physical network)."""
import functools
import requests
import pandas as pd
from string_score import uniprot_taxid, STRING, CALLER, http

INTERPRO = "https://www.ebi.ac.uk/interpro/api"
SKIP_TYPES = {"repeat", "motif", "coiled-coil", "disordered"}


def _get(url, **kw):
    r = http.get(url, timeout=60, **kw)
    r.raise_for_status()
    return r.json() if r.status_code != 204 else {}


@functools.lru_cache(maxsize=None)
def pick_pfam(acc):
    """The Pfam entry covering most residues of acc (repeats/motifs only as fallback)."""
    d = _get(f"{INTERPRO}/entry/pfam/protein/uniprot/{acc.split('-')[0]}/", params={"page_size": 200})
    cands = []
    for e in d.get("results", []):
        cov = sum(f["end"] - f["start"] + 1
                  for p in e["proteins"] for loc in p["entry_protein_locations"]
                  for f in loc["fragments"])
        cands.append((e["metadata"]["type"].lower() in SKIP_TYPES, -cov, e["metadata"]["accession"]))
    return sorted(cands)[0][2] if cands else None


@functools.lru_cache(maxsize=None)
def family_members(pfam, taxid):
    """Swiss-Prot accessions carrying `pfam` in one organism.

    InterPro first; for very large families (GPCR, WD40, SH3, ...) InterPro answers
    408 Request Timeout, so fall back to UniProt, whose Pfam xrefs come from InterPro.
    """
    url = f"{INTERPRO}/protein/reviewed/entry/pfam/{pfam}/taxonomy/uniprot/{taxid}/"
    params, out = {"page_size": 200}, []
    try:
        while url:
            d = _get(url, params=params)
            out += [x["metadata"]["accession"] for x in d.get("results", [])]
            url, params = d.get("next"), None
    except requests.HTTPError as ex:
        if ex.response is None or ex.response.status_code != 408:
            raise
        r = http.get("https://rest.uniprot.org/uniprotkb/stream",
                     params={"query": f"xref:pfam-{pfam} AND organism_id:{taxid} AND reviewed:true",
                             "fields": "accession", "format": "tsv"}, timeout=300)
        r.raise_for_status()
        out = r.text.split()[1:]
    return tuple(sorted(out))


@functools.lru_cache(maxsize=None)
def map_to_string(accs, taxid):
    """{uniprot_acc: string_id} for one organism, in one POST."""
    r = http.post(f"{STRING}/get_string_ids",
                      data={"identifiers": "\r".join(accs), "species": taxid, "limit": 1,
                            "echo_query": 1, "caller_identity": CALLER}, timeout=120)
    r.raise_for_status()
    return {h["queryItem"]: h["stringId"] for h in r.json()}


@functools.lru_cache(maxsize=None)
def string_network(string_ids, taxid, network_type):
    """All STRING edges among a frozenset of STRING ids (cached: family pairs repeat)."""
    r = http.post(f"{STRING}/network",
                      data={"identifiers": "\r".join(sorted(string_ids)), "species": taxid,
                            "required_score": 0, "network_type": network_type,
                            "caller_identity": CALLER}, timeout=300)
    r.raise_for_status()
    return tuple(r.json())


# UniProt organism ids that STRING files under another id ("unknown organism" otherwise)
STRING_TAXID = {559292: 4932}          # S. cerevisiae S288c -> S. cerevisiae


def family_hits(acc_a, acc_b, min_escore=0.4, network_type="physical"):
    """Edges between members of fam(A) and fam(B) with experiments score >= min_escore.

    Same-organism pairs only: STRING has no inter-species edges, so pairs whose proteins
    come from different organisms are returned with status "cross-species" and not searched.
    One STRING network call per family pair: the union of both families is sent and only
    cross-family edges are kept.

    Returns (summary dict, DataFrame of hits).
    """
    acc_a, acc_b = acc_a.split("-")[0], acc_b.split("-")[0]
    ta, tb = uniprot_taxid(acc_a), uniprot_taxid(acc_b)
    summ = {"A": acc_a, "B": acc_b, "taxid": ta, "taxid_B": tb}
    if ta != tb:
        summ["status"] = "cross-species"
        return summ, pd.DataFrame()
    pa, pb = pick_pfam(acc_a), pick_pfam(acc_b)
    summ.update(pfam_A=pa, pfam_B=pb)
    if not pa or not pb:
        summ["status"] = "no-pfam"
        return summ, pd.DataFrame()

    fa, fb = set(family_members(pa, ta)) | {acc_a}, set(family_members(pb, ta)) | {acc_b}
    st = STRING_TAXID.get(ta, ta)
    try:
        m = map_to_string(tuple(sorted(fa | fb)), st)
    except requests.HTTPError as ex:
        if ex.response is not None and ex.response.status_code in (400, 404):
            summ["status"] = "not-in-string"
            return summ, pd.DataFrame()
        raise
    ids_a = {m[x]: x for x in fa if x in m}
    ids_b = {m[x]: x for x in fb if x in m}
    summ.update(n_A=len(fa), n_B=len(fb), mapped_A=len(ids_a), mapped_B=len(ids_b),
                same_family=pa == pb)

    rows, seen = [], set()
    if ids_a and ids_b:
        for e in string_network(frozenset(ids_a) | frozenset(ids_b), st, network_type):
            if e["escore"] < min_escore:
                continue
            for s1, s2 in ((e["stringId_A"], e["stringId_B"]), (e["stringId_B"], e["stringId_A"])):
                if s1 in ids_a and s2 in ids_b and s1 != s2:
                    x, y = ids_a[s1], ids_b[s2]
                    if (x, y) in seen or (y, x) in seen:
                        continue
                    seen.add((x, y))
                    rows.append({"member_A": x, "gene_A": e["preferredName_A" if s1 == e["stringId_A"] else "preferredName_B"],
                                 "member_B": y, "gene_B": e["preferredName_B" if s2 == e["stringId_B"] else "preferredName_A"],
                                 "escore": e["escore"], "score": e["score"],
                                 "kind": ("A-B itself" if (x, y) == (acc_a, acc_b) else
                                          "A binds fam(B)" if x == acc_a else
                                          "B binds fam(A)" if y == acc_b else "other members")})
    hits = pd.DataFrame(rows).sort_values("escore", ascending=False) if rows else pd.DataFrame()
    summ["status"] = "ok"
    summ["n_hits"] = len(hits)
    if len(hits):
        for k in ("A-B itself", "A binds fam(B)", "B binds fam(A)", "other members"):
            summ[k] = int((hits["kind"] == k).sum())
    return summ, hits


if __name__ == "__main__":
    s, h = family_hits("P06493", "P24864")   # CDK1_HUMAN vs CCNE1_HUMAN, human

    print(s)
    print(h)
