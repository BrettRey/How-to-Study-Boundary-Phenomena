#!/usr/bin/env python3
"""Render the clarified Figure 3 from saved results, without rerunning analysis.

The original analysis script and results remain unchanged. Outputs and provenance
are written separately to data/revision/review_corrections/.
"""
from pathlib import Path
import hashlib
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/revision/results/all_specifications.csv"
OUT = ROOT / "data/revision/review_corrections"
TARGETS = ["each_other", "one_another"]
ORDER = ["none", "morph", "synt", "sem", "phon"]
MEASURES = [
    ("jaccard", "Jaccard", "#2E5077", "o"),
    ("dice", "Dice", "#C74F41", "s"),
    ("hamming", "Hamming", "#3D825D", "^"),
    ("idf_jaccard", "IDF-weighted Jaccard", "#946697", "D"),
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    style = ROOT.parents[2] / ".house-style"
    if (style / "plot_style.py").exists():
        sys.path.insert(0, str(style))
        from plot_style import setup
        setup(font_size=11, tick_size=10)
        style_source = "workspace plot_style.py"
        style_hash = sha(style / "plot_style.py")
    else:
        plt.rcParams.update({"font.family": "serif", "axes.spines.top": False,
                             "axes.spines.right": False, "legend.frameon": False})
        style_source, style_hash = "Matplotlib serif fallback", None

    result = pd.read_csv(SOURCE)
    base = result.query("pool == 'full' and holdout == 'item' and reference == 'fixed_total'")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.6), sharey=True)
    plotted = []
    for ax, item in zip(axes, TARGETS):
        for k, (metric, label, color, marker) in enumerate(MEASURES):
            rows = base.loc[base.lemma.eq(item) & base.metric.eq(metric)].set_index("drop").loc[ORDER]
            ax.scatter(np.arange(5) + (k - 1.5) * 0.12, rows.delta_percentile,
                       c=color, marker=marker, label=label, s=35)
            plotted.append(rows.reset_index()[["lemma", "pool", "holdout", "reference",
                                                "drop", "metric", "delta_percentile"]])
        ax.axhline(.5, c="0.6", ls="--", lw=.7)
        ax.set(xticks=np.arange(5),
               xticklabels=["All", "No morph.", "No syntax", "No sem.", "No phon."],
               ylim=(-.04, 1.04), title=item.replace("_", " "))
        ax.tick_params(axis="x", rotation=30)
    axes[0].set_ylabel("Reference percentile of Δ")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=10)
    fig.tight_layout(rect=(0, .13, 1, 1))
    for extension in ["pdf", "png"]:
        fig.savefig(OUT / f"specification_percentiles.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    points = pd.concat(plotted, ignore_index=True)
    points.to_csv(OUT / "plotted_values.csv", index=False)
    outputs = ["specification_percentiles.pdf", "specification_percentiles.png", "plotted_values.csv"]
    metadata = {
        "purpose": "Figure 3 legend clarification after independent review",
        "script_sha256": sha(__file__),
        "inputs": {str(SOURCE.relative_to(ROOT)): sha(SOURCE)},
        "settings": {"pool": "full", "holdout": "item", "reference": "fixed_total"},
        "metric_labels": {metric: label for metric, label, _, _ in MEASURES},
        "plotted_rows": len(points),
        "new_random_draws": 0,
        "analysis_rerun": False,
        "plot_style": style_source,
        "plot_style_sha256": style_hash,
        "versions": {"matplotlib": matplotlib.__version__, "numpy": np.__version__, "pandas": pd.__version__},
        "outputs": {name: sha(OUT / name) for name in outputs},
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Rendered Figure 3 from {len(points)} saved values: {OUT}")


if __name__ == "__main__":
    main()
