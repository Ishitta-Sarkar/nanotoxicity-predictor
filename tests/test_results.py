"""Checks that the reported numbers are the numbers produced by this code."""

import csv
import json
from pathlib import Path

import numpy as np
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

from nanotoxicity.analysis import run_analysis

ROOT = Path(__file__).resolve().parents[1]


def _metrics(y_true: np.ndarray, y_hat: np.ndarray) -> dict[str, float]:
    return {
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_hat))),
        "mae": float(mean_absolute_error(y_true, y_hat)),
        "r2": float(r2_score(y_true, y_hat)),
    }


def test_metrics_file_matches_the_prediction_table():
    payload = json.loads((ROOT / "results" / "metrics.json").read_text(encoding="utf-8"))
    rows: dict[str, list[dict[str, str]]] = {}
    with (ROOT / "results" / "loocv_predictions.csv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.setdefault(row["dataset"], []).append(row)

    by_key = {dataset["key"]: dataset for dataset in payload["datasets"]}
    assert set(by_key) == set(rows)
    for key, dataset in by_key.items():
        observed = np.array([float(row["y_observed"]) for row in rows[key]])
        predicted = np.array([float(row["y_loocv"]) for row in rows[key]])
        baseline = np.array([float(row["y_baseline"]) for row in rows[key]])
        assert len(observed) == dataset["n"]
        recomputed = {
            "model": _metrics(observed, predicted),
            "baseline": _metrics(observed, baseline),
        }
        for side in ("model", "baseline"):
            for metric in ("rmse", "mae", "r2"):
                assert abs(recomputed[side][metric] - dataset[side][metric]) < 1e-5


def test_refit_matches_the_reported_metrics(tmp_path: Path):
    payload = run_analysis(tmp_path / "cache", tmp_path / "results")
    committed = json.loads((ROOT / "results" / "metrics.json").read_text(encoding="utf-8"))
    fresh = {dataset["key"]: dataset for dataset in payload["datasets"]}
    saved = {dataset["key"]: dataset for dataset in committed["datasets"]}
    assert set(fresh) == {"gajewicz_hacat", "puzyn_ecoli"}

    gajewicz = fresh["gajewicz_hacat"]
    puzyn = fresh["puzyn_ecoli"]
    assert gajewicz["n"] == 18
    assert gajewicz["excluded_missing_endpoint"] == ["P_1", "P_2", "P_3", "P_4", "P_5"]
    assert gajewicz["descriptor_ids"] == ["dH-fc", "Xc"]
    assert puzyn["n"] == 17
    assert puzyn["excluded_missing_endpoint"] == []
    assert puzyn["descriptor_ids"] == ["dH_Me"]

    for key in fresh:
        for side in ("model", "baseline"):
            for metric in ("rmse", "mae", "r2"):
                assert abs(fresh[key][side][metric] - saved[key][side][metric]) < 1e-6
