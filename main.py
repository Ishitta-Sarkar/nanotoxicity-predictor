from src.nanotoxicity_tools import NanotoxicityPredictor


predictor = NanotoxicityPredictor(
    "data/nanomaterial_dataset.csv"
)

print("Nanotoxicity Predictor")
print("=" * 40)

summary = predictor.summarize_dataset()

print("\nDataset Summary")
print("-" * 40)
print("Total Records:", summary["total_records"])

print("\nToxicity Distribution")
for level, count in summary["toxicity_distribution"].items():
    print(f"{level}: {count}")

print("\nHighest Risk Nanomaterials")
print("-" * 40)

ranked_materials = predictor.rank_by_toxicity_risk()

for material in ranked_materials:
    print(
        material["nanomaterial"],
        "| Cell Viability:",
        material["cell_viability_percent"]
    )

print("\nHigh Toxicity Materials")
print("-" * 40)

high_toxicity = predictor.filter_high_toxicity()

for material in high_toxicity:
    print(material["nanomaterial"])

print("\nRisk Summary")
print("-" * 40)

print(predictor.generate_risk_summary())
