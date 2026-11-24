"""STRING score for a pair of UniProt accessions (experiments channel by default)."""
import functools
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# one session for all calls: retries 429/5xx with exponential backoff (2, 4, 8, 16, 32 s)
http = requests.Session()
http.mount("https://", HTTPAdapter(max_retries=Retry(
    total=5, backoff_factor=2, status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=None, respect_retry_after_header=True)))

STRING = "https://string-db.org/api/json"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
CALLER = "python-client"


@functools.lru_cache(maxsize=None)
def uniprot_taxid(acc):
    """NCBI taxon id of a UniProt accession (isoform suffix dropped)."""
    r = http.get(f"{UNIPROT}/{acc.split('-')[0]}.json",
                     params={"fields": "organism_id"}, timeout=30)
    r.raise_for_status()
    return r.json()["organism"]["taxonId"]


@functools.lru_cache(maxsize=None)
def string_id(acc, taxid):
    """Map a UniProt accession to a STRING id within one species, or None."""
    r = http.get(f"{STRING}/get_string_ids",
                     params={"identifiers": acc.split("-")[0], "species": taxid,
                             "limit": 1, "caller_identity": CALLER}, timeout=30)
    r.raise_for_status()
    hits = r.json()
    return hits[0]["stringId"] if hits else None


def string_score(acc_a, acc_b, channel="escore", network_type="physical"):
    """Return the STRING score (0-1) between two UniProt accessions.

    channel: 'escore' (experiments, default), 'dscore' (databases),
             'tscore' (text mining), 'score' (combined), ...
    Returns 0.0 if both map but STRING has no edge,
            None if the pair can't be scored (different species or unmapped id).
    """
    ta, tb = uniprot_taxid(acc_a), uniprot_taxid(acc_b)
    if ta != tb:                      # STRING has no cross-species edges
        return None
    sa, sb = string_id(acc_a, ta), string_id(acc_b, tb)
    if not sa or not sb:
        return None
    r = http.get(f"{STRING}/network",
                     params={"identifiers": f"{sa}%0d{sb}", "species": ta,
                             "required_score": 0, "network_type": network_type,
                             "caller_identity": CALLER}, timeout=30)
    r.raise_for_status()
    for e in r.json():
        if {e["stringId_A"], e["stringId_B"]} == {sa, sb}:
            return e[channel]
    return 0.0


if __name__ == "__main__":
    for a, b in [("P24941", "P24864"),    # CDK2-CCNE1, known binders
                 ("P62879", "Q92730"),    # GNB2-RND1, Negatome
                 ("P43407", "P21981")]:   # Sdc2-Tgm2 (mouse), Negatome
        print(a, b, string_score(a, b), string_score(a, b, "score"))
