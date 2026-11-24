"""Run the download, the leave-one-out fits, and the result files."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from nanotoxicity.datasets import ARCHIVES
from nanotoxicity.model import leave_one_out
from nanotoxicity.plot import save_figure
from nanotoxicity.qdb import OxideTable, download_archive, read_archive


def run_analysis(cache_dir: Path, results_dir: Path) -> dict:
    tables = []
    fits = []
    for spec in ARCHIVES:
        archive_path = download_archive(spec, cache_dir)
        table = read_archive(archive_path, spec)
        tables.append(table)
        fits.append(leave_one_out(table.x, table.y))

    results_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "validation": "leave-one-out",
        "model": "ordinary least squares on descriptors standardized within each training fold",
        "baseline": "leave-one-out mean of the endpoint",
        "datasets": [_dataset_record(table, fit) for table, fit in zip(tables, fits)],
    }
    (results_dir / "metrics.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    _write_predictions(results_dir / "loocv_predictions.csv", tables, fits)
    save_figure(tables, fits, results_dir / "loocv_observed_vs_predicted.png")
    return payload


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Leave-one-out reanalysis of two published metal-oxide nano-QSAR tables."
    )
    parser.add_argument("--cache", type=Path, default=Path("data/cache"))
    parser.add_argument("--results", type=Path, default=Path("results"))
    args = parser.parse_args(argv)
    payload = run_analysis(args.cache, args.results)
    print("Leave-one-out results")
    print("Descriptors are the columns deposited as descriptors. The endpoint is not an input.")
    print()
    for dataset in payload["datasets"]:
        print(f"{dataset['citation']} — {dataset['system']} (n = {dataset['n']})")
        print(f"  endpoint: {dataset['endpoint']}")
        print(f"  descriptors: {', '.join(dataset['descriptor_ids'])}")
        if dataset["excluded_missing_endpoint"]:
            print(
                "  excluded, no measured endpoint: "
                + ", ".join(dataset["excluded_missing_endpoint"])
            )
        model = dataset["model"]
        baseline = dataset["baseline"]
        print(
            f"  model     RMSE {model['rmse']:.3f}   MAE {model['mae']:.3f}   R2 {model['r2']:.3f}"
        )
        print(
            "  baseline  "
            f"RMSE {baseline['rmse']:.3f}   MAE {baseline['mae']:.3f}   R2 {baseline['r2']:.3f}"
        )
        print(f"  wrote {args.results}")
        print()


def _dataset_record(table: OxideTable, fit) -> dict:
    return {
        "key": table.key,
        "citation": table.citation,
        "system": table.system,
        "n": len(table.y),
        "endpoint": table.endpoint_label,
        "property_id": table.property_id,
        "descriptor_ids": list(table.descriptor_ids),
        "descriptor_names": list(table.descriptor_names),
        "descriptor_units": list(table.descriptor_units),
        "excluded_missing_endpoint": list(table.excluded_missing_endpoint),
        "model": {
            "rmse": fit.model.rmse,
            "mae": fit.model.mae,
            "r2": fit.model.r2,
        },
        "baseline": {
            "rmse": fit.baseline.rmse,
            "mae": fit.baseline.mae,
            "r2": fit.baseline.r2,
        },
        "full_sample_standardized_coefficients": {
            "intercept": fit.intercept,
            "per_descriptor": {
                descriptor_id: coefficient
                for descriptor_id, coefficient in zip(table.descriptor_ids, fit.coefficients)
            },
            "note": (
                "Fit on every labelled row after standardization. "
                "These coefficients are descriptive. They are not cross-validated."
            ),
        },
    }


def _write_predictions(path: Path, tables: list[OxideTable], fits: list) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "dataset",
                "compound_id",
                "name",
                "formula",
                "y_observed",
                "y_loocv",
                "y_baseline",
                "descriptors",
            ]
        )
        for table, fit in zip(tables, fits):
            for index, compound_id in enumerate(table.compound_ids):
                descriptor_text = ";".join(
                    f"{descriptor_id}={table.x[index, column]:.6g}"
                    for column, descriptor_id in enumerate(table.descriptor_ids)
                )
                writer.writerow(
                    [
                        table.key,
                        compound_id,
                        table.names[index],
                        table.formulas[index],
                        f"{table.y[index]:.6f}",
                        f"{fit.y_pred[index]:.6f}",
                        f"{fit.y_baseline[index]:.6f}",
                        descriptor_text,
                    ]
                )
