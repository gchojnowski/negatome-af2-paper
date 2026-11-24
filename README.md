# Can AlphaFold2 say no? Interface confidence on curated non-interacting proteins

Data, analysis code and figure scripts accompanying the preprint:

> Chojnowski, G. (2026). *Can AlphaFold2 say no? Interface confidence on curated non-interacting proteins.*

## Overview

AlphaFold2 complex predictions were generated for 1,000 heteromeric protein pairs annotated as non-interacting in the Negatome 2.0 database and scored with ipTM and ipSAE. The literature evidence behind each annotation was re-assessed with LLM-assisted literature review, and evolutionary relationships between proteins were defined using Pfam families. The analysis shows that high-confidence predictions for non-interacting pairs are enriched among proteins whose homologues are known to interact.

## Repository contents

```
manuscript/   Preprint PDF
data/         Pair lists, prediction scores, evidence classes, Pfam assignments
analysis/     Notebooks and scripts for the statistical analyses
figures/      Scripts that regenerate the figures in the manuscript
```

## Citation

If you use this data or code, please cite the preprint above. Citation metadata is provided in `CITATION.cff`; the archived version of this repository is available on Zenodo: DOI:

## License

- Code: MIT License
- Manuscript, data and figures: Creative Commons Attribution 4.0 International (CC BY 4.0)

See `LICENSE` for details.

## Contact

Grzegorz Chojnowski, Independent researcher, Hamburg, Germany
ORCID: 0000-0002-3796-8352
