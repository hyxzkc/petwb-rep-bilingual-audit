# PETWB-REP bilingual report audit

Reproduction scripts and locked analysis outputs for an audit of clinical semantic drift between paired Chinese and English whole-body FDG-PET/CT reports in PETWB-REP (490 cases). The study codes organ-level uptake states, separates rule-screening signals from extraction artifacts, and reports stratum-weighted PPV together with inter-annotator agreement.

## Data

Raw reports are not included here. Download them from Zenodo and check MD5 against `data_version_lock.md` before running anything:

- English reports and metadata: DOI 10.5281/zenodo.18670487 (`report_en.csv`, `meta_data.csv`)
- Chinese reports: DOI 10.5281/zenodo.19649618 (`report_zh.xlsx`, converted to UTF-8 CSV)

The v2 Chinese CSV is corrupted in the released file and is kept only as a data-quality artifact. The analysis used v3 Chinese with v2 English and v2 metadata; the pairing was verified row by row on Image ID.

## Reproduce the PPV

`reproduce_ppv_bootstrap.py` uses only the Python standard library. Two modes are available:

```
python reproduce_ppv_bootstrap.py --mode structured
python reproduce_ppv_bootstrap.py --mode expert
```

Both apply the six prespecified strata and their population weights to the 120 candidates drawn from the 2,012 over-normalization candidates, then resample within strata 20,000 times with seed 20261120. The structured mode uses the earlier structured-adjudication labels; the expert mode uses the blinded review by a nuclear medicine physician. Each run prints the weighted PPV, per-stratum counts, and a bootstrap 95% interval. The 60 rule-match controls belong to the agreement set only and are not pooled into PPV.

`--seed` and `--replicates` can be overridden, but the defaults are the locked settings.

## Files

- `reproduce_ppv_bootstrap.py`: PPV and interval reproduction
- `data_version_lock.md`: source versions, DOIs, MD5, and pairing checks
- `organ_state_extractor.py`: rule-based organ and uptake-state extraction
- `organ_state_extractor_pre_fix.py`: pre-fix extractor, kept so the evidence-length change remains auditable
- `overnorm_validation_sample120_ppv.csv`: the 120 stratified candidates
- `overnorm120_adjudicated.json`, `adjudication_all338.json`: structured adjudication records
- `expert_gold_standard.json`: blinded expert labels for all 180 units
- `overnorm_ppv_result.json`, `expert_ppv_result.json`, `full_stats_adjudicated.json`, `kappa_improved_result.json`, `expert_agreement_report.md`: locked results

Do not rerun the extractor over the locked outputs; write to a separate directory and version any new results.

## Data use and citation

Use the PETWB-REP releases under their published terms. No re-identification and no clinical use. Cite the dataset paper (Scientific Data, DOI 10.1038/s41597-026-07058-w) and the Zenodo releases above.
