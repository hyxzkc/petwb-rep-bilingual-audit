#!/usr/bin/env python3
"""Reproduce stratified over-normalization PPV estimates.

Uses the locked 120-row sample and six sampling strata. The structured mode
reads the adjudicated labels. The expert mode reads E-0731 labels from the
expert gold-standard file. The script uses only the Python standard library.
"""
import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SAMPLE = ROOT / "overnorm_validation_sample120_ppv.csv"
RESULT = ROOT / "overnorm_ppv_result.json"
EXPERT_RESULT = ROOT / "expert_ppv_result.json"
STRUCTURED = ROOT / "overnorm120_adjudicated.json"
EXPERT = ROOT / "expert_gold_standard.json"


def load_labels(mode):
    if mode == "structured":
        rows = json.loads(STRUCTURED.read_text(encoding="utf-8"))
        return [row["final_label"] for row in rows]
    rows = json.loads(EXPERT.read_text(encoding="utf-8"))[:120]
    return [row["expert_label"] for row in rows]


def load_sample():
    with SAMPLE.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def calculate(sample, labels, seed=20261120, replicates=20000):
    if len(sample) != 120 or len(labels) != 120:
        raise ValueError("The PPV sample and labels must each contain 120 rows")
    population = {
        "not_abnormal|plain": 1499,
        "not_abnormal|quantified_or_qualified": 381,
        "not_increased|plain": 118,
        "not_increased|quantified_or_qualified": 4,
        "other|plain": 8,
        "other|quantified_or_qualified": 2,
    }
    strata = defaultdict(list)
    for index, row in enumerate(sample):
        strata[row["sampling_stratum"]].append(index)
    details = {}
    for name, indices in strata.items():
        positives = sum(labels[index] == "over_normalization" for index in indices)
        n = len(indices)
        details[name] = {
            "N": population[name],
            "n": n,
            "over_norm": positives,
            "ppv": positives / n,
            "weight": population[name] / n,
        }
    weighted = sum(item["N"] * item["ppv"] for item in details.values()) / 2012
    rng = random.Random(seed)
    boot = []
    for _ in range(replicates):
        numerator = 0.0
        for name, indices in strata.items():
            values = [labels[index] == "over_normalization" for index in indices]
            resampled = [values[rng.randrange(len(values))] for _ in values]
            numerator += population[name] * (sum(resampled) / len(resampled))
        boot.append(numerator / 2012)
    boot.sort()
    low = boot[int(replicates * 0.025)]
    high = boot[int(replicates * 0.975) - 1]
    return details, weighted, (low, high)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["structured", "expert"], default="expert")
    parser.add_argument("--seed", type=int, default=20261120)
    parser.add_argument("--replicates", type=int, default=20000)
    args = parser.parse_args()
    sample = load_sample()
    labels = load_labels(args.mode)
    details, estimate, ci = calculate(sample, labels, args.seed, args.replicates)
    stored = json.loads((RESULT if args.mode == "structured" else EXPERT_RESULT).read_text(encoding="utf-8"))
    print(json.dumps({
        "mode": args.mode,
        "estimate": estimate,
        "percent": f"{estimate * 100:.2f}%",
        "bootstrap_95ci": ci,
        "ci_percent": f"{ci[0] * 100:.2f}%–{ci[1] * 100:.2f}%",
        "seed": args.seed,
        "replicates": args.replicates,
        "stored_expert_estimate": stored.get("ppv", {}).get("estimate"),
        "stored_expert_ci": stored.get("ppv", {}).get("bootstrap_95ci"),
        "strata": details,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
