#!/usr/bin/env python3
"""Approved post-review block weighting and summaries of existing results.

Writes only data/revision/followup; does not rerun or overwrite the original analysis.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import revision_analysis as original

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/revision/followup"
SAVED = ROOT / "data/revision/results"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def block_weights(blocks):
    blocks = np.asarray(blocks)
    return np.array([1.0 / np.count_nonzero(blocks == b) for b in blocks])


def weighted_jaccard(target, anchors, weights):
    target = np.asarray(target, dtype=bool)
    anchors = np.asarray(anchors, dtype=bool)
    weights = np.asarray(weights, dtype=float)
    numerator = ((anchors != target) * weights).sum(axis=1)
    denominator = ((anchors | target) * weights).sum(axis=1)
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)


def signed_margin(frame):
    return np.where(frame.role.eq("pronoun_anchor"), -frame.delta,
                    np.where(frame.role.eq("compound_anchor"), frame.delta, np.nan))


def summaries(frame, keys):
    rows = []
    for key, group in frame.groupby(keys, sort=True):
        if not isinstance(key, tuple):
            key = (key,)
        controls = group[group.eligible_control]
        pron = controls[controls.role.eq("pronoun_anchor")].delta
        comp = controls[controls.role.eq("compound_anchor")].delta
        row = dict(zip(keys, key))
        row.update(n_controls=len(controls), controls_correct=int((signed_margin(controls) > 0).sum()),
                   p_min=pron.min(), p_max=pron.max(), c_min=comp.min(), c_max=comp.max())
        for item in original.TARGETS:
            target = group[group.lemma.eq(item)].iloc[0]
            row[item + "_delta"] = target.delta
            row[item + "_within_pronoun_range"] = bool(pron.min() <= target.delta <= pron.max())
            row[item + "_within_compound_range"] = bool(comp.min() <= target.delta <= comp.max())
            row[item + "_between_control_ranges"] = bool(pron.max() < target.delta < comp.min())
        rows.append(row)
    return pd.DataFrame(rows)


def new_weighting(x, manifest, blocks):
    records, weights_out = [], []
    for drop in ["none", "morph"]:
        keep = np.ones(len(blocks), dtype=bool) if drop == "none" else blocks.ne("morph").to_numpy()
        xx, bb = x.loc[:, keep], blocks.loc[keep]
        weights = block_weights(bb)
        for feature, block, weight in zip(xx.columns, bb, weights):
            weights_out.append(dict(drop=drop, feature=feature, block=block, weight=weight))
        for holdout in ["item", "family"]:
            for item in manifest.index[manifest.role.ne("context")]:
                mask = original.comparator_mask(manifest, item, "full", holdout)
                anchors = xx.loc[mask]
                pron = manifest.loc[mask, "role"].eq("pronoun_anchor").to_numpy()
                distances = weighted_jaccard(xx.loc[item], anchors, weights)
                dp, dc = distances[pron].mean(), distances[~pron].mean()
                role = manifest.loc[item, "role"]
                records.append(dict(lemma=item, role=role, family=manifest.loc[item, "family"],
                    eligible_control=role != "target", pool="full", holdout=holdout, drop=drop,
                    weighting="equal_block_total", features=xx.shape[1],
                    n_pronoun=int(pron.sum()), n_compound=int((~pron).sum()),
                    d_pronoun=dp, d_compound=dc, delta=dp-dc))
    pd.DataFrame(weights_out).to_csv(OUT / "feature_weights.csv", index=False)
    result = pd.DataFrame(records)
    result["signed_control_margin"] = signed_margin(result)
    result.to_csv(OUT / "block_weight_item_results.csv", index=False)
    return result


def existing_reporting(saved):
    keys = ["pool", "holdout", "drop", "metric", "reference"]
    grid = summaries(saved, keys)
    target_ranks = saved[saved.role.eq("target")].pivot(index=keys, columns="lemma", values="delta_percentile")
    target_ranks.columns = [item + "_reference_percentile" for item in target_ranks.columns]
    grid = grid.merge(target_ranks.reset_index(), on=keys, validate="one_to_one")
    assert len(grid) == 160
    grid.to_csv(OUT / "complete_specification_grid.csv", index=False)
    # Observed distances do not depend on the placement benchmark: avoid duplicating them.
    controls = saved[saved.reference.eq("fixed_total") & saved.eligible_control].copy()
    controls["signed_control_margin"] = signed_margin(controls)
    cols = ["lemma", "role", "family", "pool", "holdout", "drop", "metric",
            "d_pronoun", "d_compound", "delta", "signed_control_margin"]
    controls[cols].to_csv(OUT / "original_control_margins.csv", index=False)
    group = controls.groupby(["pool", "holdout", "drop", "metric", "role"], sort=True)
    margin_summary = group.signed_control_margin.agg(count="size", minimum="min", median="median", maximum="max").reset_index()
    margin_summary.to_csv(OUT / "original_control_margin_summary.csv", index=False)

    draws = pd.read_csv(SAVED / "target_reference_draws.csv.gz")
    rows = []
    for (item, reference), group in draws.groupby(["lemma", "reference"], sort=True):
        obs = saved.query("lemma==@item and reference==@reference and pool=='full' and holdout=='item' and drop=='none' and metric=='jaccard'").iloc[0]
        for quantity in ["d_pronoun", "d_compound", "delta"]:
            values = group[quantity].to_numpy()
            summary = original.rank_summary(float(obs[quantity]), values)
            below = int((values < obs[quantity] - original.TOL).sum())
            tied = int((np.abs(values - obs[quantity]) <= original.TOL).sum())
            for name in ["percentile", "percentile_mcse", "lower_tail"]:
                assert np.isclose(summary[name], obs[quantity+"_"+name], atol=1e-12, rtol=0)
            rows.append(dict(lemma=item, reference=reference, quantity=quantity, observed=obs[quantity],
                draws=len(values), strictly_below=below, tied=tied, lower_tail_hits=below+tied,
                **{name:summary[name] for name in ["percentile", "percentile_mcse", "lower_tail"]},
                zero_hit_probability_upper_95=(-np.expm1(np.log(.05)/len(values)) if below+tied == 0 else np.nan)))
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "target_monte_carlo_details.csv", index=False)
    return result


def tables(comparison, summary, mc):
    rows = []
    for drop, label in [("none", "All features"), ("morph", "No morphology")]:
        for weighting, name in [("original_equal_feature", "Equal feature"), ("equal_block_total", "Equal block total")]:
            s = summary.query("drop==@drop and weighting==@weighting").set_index("holdout")
            rows.append(label+" & "+name+" & "+" & ".join(f"{s.loc['item',item+'_delta']:.2f}" for item in original.TARGETS)+
                        " & "+" & ".join(f"{int(s.loc[h,'controls_correct'])}/{int(s.loc[h,'n_controls'])}" for h in ["item", "family"])+r" \\")
    (OUT / "block_weight_table.tex").write_text("\n".join(rows)+"\n"+r"\bottomrule"+"\n")
    macros = {}
    for item, prefix in [("each_other", "Each"), ("one_another", "One")]:
        z = comparison.query("weighting=='equal_block_total' and drop=='none' and holdout=='item' and lemma==@item").iloc[0]
        for quantity, suffix in [("d_pronoun", "DP"), ("d_compound", "DC"), ("delta", "Delta")]:
            macros[prefix+suffix] = f"{z[quantity]:.2f}"
        z = mc.query("lemma==@item and reference=='fixed_total' and quantity=='delta'").iloc[0]
        macros[prefix+"MCSE"] = f"{z.percentile_mcse:.2g}"
    macros["ZeroHitUpper"] = f"{mc.zero_hit_probability_upper_95.dropna().iloc[0]:.2g}"
    (OUT / "numbers.tex").write_text("% Generated by code/revision_followup.py; do not hand-edit.\n"+
        "\n".join(r"\newcommand{\Follow"+name+"}{"+value+"}" for name,value in macros.items())+"\n")


def main():
    protocol_path = ROOT / "data/revision/followup_protocol.json"
    protocol = json.loads(protocol_path.read_text())
    for name, digest in protocol["inputs"].items():
        if sha(ROOT / name) != digest:
            raise ValueError("Input differs from approved follow-up protocol: " + name)
    OUT.mkdir(exist_ok=True)
    _, x, manifest, blocks, _ = original.load_design()
    saved = pd.read_csv(SAVED / "all_specifications.csv")
    result = new_weighting(x, manifest, blocks)
    assert len(result) == 236 and result.groupby(["drop", "holdout"]).size().eq(59).all()
    baseline = saved.query("metric=='jaccard' and pool=='full' and drop in ['none','morph'] and reference=='fixed_total'").copy()
    baseline["weighting"] = "original_equal_feature"
    baseline["signed_control_margin"] = signed_margin(baseline)
    comparison = pd.concat([baseline[result.columns], result], ignore_index=True)
    comparison.to_csv(OUT / "weighting_comparison.csv", index=False)
    summary = summaries(comparison, ["drop", "holdout", "weighting"])
    summary.to_csv(OUT / "block_weight_summary.csv", index=False)
    mc = existing_reporting(saved)
    tables(comparison, summary, mc)
    metadata = dict(protocol_sha256=sha(protocol_path), script_sha256=sha(__file__),
                    original_analysis_script_sha256=sha(ROOT / "code/revision_analysis.py"),
                    inputs=protocol["inputs"], numpy=np.__version__, pandas=pd.__version__,
                    new_weighted_rows=len(result), comparison_rows=len(comparison), new_random_draws=0,
                    outputs={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name != "metadata.json"})
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
