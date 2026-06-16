import csv


class NanotoxicityPredictor:
    def __init__(self, file_path):
        self.file_path = file_path
        self.records = self.load_data()

    def load_data(self):
        records = []

        with open(self.file_path, "r", encoding="utf-8") as file:
            reader = csv.DictReader(file)
            for row in reader:
                records.append(row)

        return records

    def summarize_dataset(self):
        total_records = len(self.records)
        toxicity_counts = {}

        for record in self.records:
            level = record["toxicity_level"]
            toxicity_counts[level] = toxicity_counts.get(level, 0) + 1

        return {
            "total_records": total_records,
            "toxicity_distribution": toxicity_counts
        }

    def rank_by_toxicity_risk(self):
        ranked = sorted(
            self.records,
            key=lambda record: float(record["cell_viability_percent"])
        )

        return ranked

    def filter_high_toxicity(self):
        return [
            record for record in self.records
            if record["toxicity_level"].lower() == "high"
        ]

    def estimate_risk_score(self, record):
        size = float(record["particle_size_nm"])
        zeta = abs(float(record["zeta_potential_mv"]))
        surface_area = float(record["surface_area_m2g"])
        concentration = float(record["concentration_ugml"])
        viability = float(record["cell_viability_percent"])

        risk_score = (
            (100 - viability) * 0.4
            + concentration * 0.25
            + surface_area * 0.2
            + zeta * 0.1
            + (100 / size) * 0.05
        )

        return round(risk_score, 2)

    def generate_risk_summary(self):
        summary = "Nanomaterial Toxicity Risk Summary\n"
        summary += "=" * 45 + "\n\n"

        for record in self.records:
            score = self.estimate_risk_score(record)
            summary += (
                f"{record['nanomaterial']} | "
                f"Toxicity Level: {record['toxicity_level']} | "
                f"Estimated Risk Score: {score}\n"
            )

        return summary.strip()
