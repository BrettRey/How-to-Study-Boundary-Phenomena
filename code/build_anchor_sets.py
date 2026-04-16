from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "inventory_annotations.csv"
CSV_OUTPUT = ROOT / "data" / "anchor_sets.csv"
PRON_OUTPUT = ROOT / "data" / "clear_pronoun_anchors.txt"
DET_OUTPUT = ROOT / "data" / "clear_determinative_anchors.txt"
REVIEW_OUTPUT = ROOT / "data" / "anchor_review_items.txt"
EXCLUDE_OUTPUT = ROOT / "data" / "anchor_excluded_items.txt"


def anchor_bucket(row: dict[str, str]) -> str:
    recommendation = row["anchor_recommendation"]
    revised_class = row["revised_class"]
    if recommendation == "retain" and revised_class == "pronoun":
        return "clear_pronoun_anchor"
    if recommendation == "retain" and revised_class == "determinative":
        return "clear_determinative_anchor"
    if recommendation == "review":
        return "review_item"
    if recommendation == "exclude":
        return "excluded_item"
    raise ValueError(
        f"Unsupported anchor mapping for lemma={row['lemma']}: "
        f"recommendation={recommendation}, revised_class={revised_class}"
    )


def write_manifest(path: Path, items: list[str]) -> None:
    path.write_text("\n".join(items) + ("\n" if items else ""))


def main() -> None:
    with INPUT.open(newline="") as f:
        rows = list(csv.DictReader(f))

    for row in rows:
        row["anchor_bucket"] = anchor_bucket(row)

    with CSV_OUTPUT.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "lemma",
                "revised_class",
                "anchor_bucket",
                "category_source",
                "personhood_profile",
                "notes",
            ],
        )
        writer.writeheader()
        writer.writerows(
            {
                "lemma": row["lemma"],
                "revised_class": row["revised_class"],
                "anchor_bucket": row["anchor_bucket"],
                "category_source": row["category_source"],
                "personhood_profile": row["personhood_profile"],
                "notes": row["notes"],
            }
            for row in rows
        )

    clear_pronouns = [
        row["lemma"] for row in rows if row["anchor_bucket"] == "clear_pronoun_anchor"
    ]
    clear_determinatives = [
        row["lemma"]
        for row in rows
        if row["anchor_bucket"] == "clear_determinative_anchor"
    ]
    review_items = [row["lemma"] for row in rows if row["anchor_bucket"] == "review_item"]
    excluded_items = [
        row["lemma"] for row in rows if row["anchor_bucket"] == "excluded_item"
    ]

    write_manifest(PRON_OUTPUT, clear_pronouns)
    write_manifest(DET_OUTPUT, clear_determinatives)
    write_manifest(REVIEW_OUTPUT, review_items)
    write_manifest(EXCLUDE_OUTPUT, excluded_items)


if __name__ == "__main__":
    main()
