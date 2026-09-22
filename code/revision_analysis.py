#!/usr/bin/env python3
"""Fixed-anchor distance analysis. Run from any directory; no input is rewritten.

The protocol and comparator manifest in data/revision define the design.
Requires numpy, pandas and matplotlib. Outputs go to data/revision/results.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TARGETS = ("each_other", "one_another")
METRICS = ("jaccard", "dice", "hamming", "idf_jaccard")
TOL = 1e-12


def load_design():
    protocol = json.loads((ROOT / "data/revision/protocol.json").read_text())
    for name, digest in protocol["inputs"].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Input changed after protocol was recorded: {name}")
    raw = pd.read_csv(ROOT / "data/matrix_clean.csv").set_index("lemma")
    manifest = pd.read_csv(ROOT / "data/revision/comparators.csv").set_index("lemma")
    assert raw.index.is_unique and manifest.index.is_unique
    assert raw.index.equals(manifest.index)
    full = raw.drop(columns="class")
    assert np.isin(full.to_numpy(), [0, 1]).all()
    retained = full.columns[full.sum(axis=0) >= 2]
    excluded = full.columns.difference(retained).tolist()
    assert excluded == ["each_other"], excluded
    blocks = pd.read_csv(ROOT / "data/feature_blocks_by_order.csv").set_index("feature")
    assert blocks.index.is_unique and set(blocks.index) == set(full.columns)
    assert set(manifest.loc[list(TARGETS), "role"]) == {"target"}
    return protocol, full[retained], manifest, blocks.loc[retained, "block"], excluded


def reference_profiles(x, blocks, mode, count, rng):
    """Uniform subset of positions, conditional on total or blockwise totals."""
    z = np.zeros((count, len(x)), dtype=float)
    groups = [np.arange(len(x))] if mode == "fixed_total" else [
        np.flatnonzero(blocks == b) for b in sorted(set(blocks))
    ]
    for cols in groups:
        k = int(x[cols].sum())
        if k == 0:
            continue
        if k == len(cols):
            z[:, cols] = 1
            continue
        # Independent continuous ranks give every k-subset equal probability.
        ranks = rng.random((count, len(cols)))
        selected = np.argpartition(ranks, k - 1, axis=1)[:, :k]
        z[np.arange(count)[:, None], cols[selected]] = 1
    assert np.all(z.sum(axis=1) == x.sum())
    return z


def distance_matrices(z, a):
    """Rows of z to rows of a; weights estimated from a alone."""
    z = np.asarray(z, dtype=float)
    a = np.asarray(a, dtype=float)
    shared = z @ a.T
    total = z.sum(axis=1)[:, None] + a.sum(axis=1)[None, :]
    union = total - shared
    yield "jaccard", np.divide(total - 2 * shared, union,
                              out=np.zeros_like(shared), where=union > 0)
    yield "dice", np.divide(total - 2 * shared, total,
                           out=np.zeros_like(shared), where=total > 0)
    yield "hamming", (total - 2 * shared) / z.shape[1]
    weights = np.log(len(a) / np.maximum(1, a.sum(axis=0)))
    shared_w = (z * weights) @ a.T
    total_w = (z @ weights)[:, None] + (a @ weights)[None, :]
    union_w = total_w - shared_w
    dw = np.divide(total_w - 2 * shared_w, union_w,
                   out=np.zeros_like(shared_w), where=union_w > TOL)
    yield "idf_jaccard", np.clip(dw, 0, 1)


def jaccard(z, a):
    return next(distance_matrices(np.atleast_2d(z), np.atleast_2d(a)))[1]


def comparator_mask(manifest, item, pool, holdout):
    mask = manifest.role.isin(["pronoun_anchor", "compound_anchor"])
    if pool == "nonlocative":
        mask &= manifest.locative_compound.eq(0)
    # A locative compound omitted by the pool remains a held-out probe, not a control.
    if holdout == "family" and manifest.loc[item, "role"] != "target":
        mask &= manifest.family.ne(manifest.loc[item, "family"])
    else:
        mask.loc[item] = False
    assert not mask.loc[list(TARGETS)].any() and not mask.loc[item]
    return mask


def rank_summary(obs, null):
    below = np.count_nonzero(null < obs - TOL)
    tied = np.count_nonzero(np.abs(null - obs) <= TOL)
    q = (below + 0.5 * tied) / len(null)
    lower = (1 + below + tied) / (len(null) + 1)
    quant = np.quantile(null, [0.025, 0.5, 0.975])
    # SE of the empirical midrank mean (including ties), not a sampling SE over words.
    score = (null < obs - TOL).astype(float) + 0.5 * (np.abs(null - obs) <= TOL)
    se = score.std(ddof=1) / np.sqrt(len(null))
    return dict(percentile=q, percentile_mcse=se, lower_tail=lower,
                null_q025=quant[0], null_median=quant[1], null_q975=quant[2])


def analyse(x, manifest, blocks, protocol, draws, out):
    results = []
    primary_draws = []
    active = manifest.index[manifest.role.ne("context")]
    for drop in protocol["ablations"]:
        keep = np.ones(len(blocks), dtype=bool) if drop == "none" else blocks.ne(drop).to_numpy()
        xx = x.loc[:, keep]
        bb = blocks[keep].to_numpy()
        print(f"Analysing drop={drop}, {xx.shape[1]} features", flush=True)
        for item in active:
            target = xx.loc[item].to_numpy(dtype=float)
            for mode in protocol["references"]:
                token = f"{protocol['seed']}:{drop}:{item}:{mode}".encode()
                seed = int.from_bytes(hashlib.sha256(token).digest()[:8], "little")
                null_x = reference_profiles(target, bb, mode, draws, np.random.default_rng(seed))
                z = np.vstack([target, null_x])
                for pool in protocol["comparator_pools"]:
                    for holdout in protocol["holdouts"]:
                        mask = comparator_mask(manifest, item, pool, holdout)
                        a = xx.loc[mask].to_numpy(dtype=float)
                        pron = manifest.loc[mask, "role"].eq("pronoun_anchor").to_numpy()
                        is_control = manifest.loc[item, "role"] in ["pronoun_anchor", "compound_anchor"]
                        if pool == "nonlocative" and manifest.loc[item, "locative_compound"]:
                            is_control = False
                        for metric, d in distance_matrices(z, a):
                            d_p, d_c = d[:, pron].mean(axis=1), d[:, ~pron].mean(axis=1)
                            stats = {"d_pronoun": d_p, "d_compound": d_c, "delta": d_p - d_c}
                            row = dict(lemma=item, role=manifest.loc[item, "role"],
                                       eligible_control=is_control, family=manifest.loc[item, "family"],
                                       pool=pool, holdout=holdout, drop=drop, metric=metric,
                                       reference=mode, features=xx.shape[1], ones=int(target.sum()),
                                       n_pronoun=int(pron.sum()), n_compound=int((~pron).sum()), draws=draws)
                            for name, values in stats.items():
                                row[name] = values[0]
                                for key, val in rank_summary(values[0], values[1:]).items():
                                    row[f"{name}_{key}"] = val
                            results.append(row)
                            if item in TARGETS and drop == "none" and pool == "full" and holdout == "item" and metric == "jaccard":
                                primary_draws.append(pd.DataFrame({"lemma": item, "reference": mode,
                                    "draw": np.arange(1, draws + 1), **{k:v[1:] for k,v in stats.items()}}))
    result = pd.DataFrame(results)
    result.to_csv(out / "all_specifications.csv", index=False)
    primary = result.query("drop == 'none' and pool == 'full' and holdout == 'item' and metric == 'jaccard' and reference == 'fixed_total'")
    primary.to_csv(out / "primary_item_results.csv", index=False)
    pd.concat(primary_draws, ignore_index=True).to_csv(out / "target_reference_draws.csv.gz", index=False,
                                                    compression={"method":"gzip", "mtime":0})
    return result, primary


def coding_sensitivity(x, manifest, out):
    anchor = manifest.role.isin(["pronoun_anchor", "compound_anchor"])
    a = x.loc[anchor].to_numpy(dtype=float)
    pron = manifest.loc[anchor, "role"].eq("pronoun_anchor").to_numpy()
    rows = []
    for item in TARGETS:
        target = x.loc[item].to_numpy(dtype=float)
        for feature_idx, feature in enumerate(x.columns):
            for operation in ["delete_feature", "reverse_polarity", "flip_target_cell"]:
                aa, tt = a.copy(), target.copy()
                if operation == "delete_feature":
                    aa, tt = np.delete(aa, feature_idx, axis=1), np.delete(tt, feature_idx)
                elif operation == "reverse_polarity":
                    aa[:, feature_idx] = 1 - aa[:, feature_idx]
                    tt[feature_idx] = 1 - tt[feature_idx]
                else:
                    tt[feature_idx] = 1 - tt[feature_idx]
                for metric, d in distance_matrices(tt[None, :], aa):
                    dp, dc = d[0, pron].mean(), d[0, ~pron].mean()
                    rows.append(dict(lemma=item, feature=feature, operation=operation, metric=metric,
                                     d_pronoun=dp, d_compound=dc, delta=dp-dc))
    pd.DataFrame(rows).to_csv(out / "coding_sensitivity.csv", index=False)


def resample_families(x, manifest, seed, out):
    rng = np.random.default_rng(seed)
    anchors = manifest.role.isin(["pronoun_anchor", "compound_anchor"])
    aa = x.loc[anchors]
    mm = manifest.loc[anchors]
    groups = {role: list(mm.loc[mm.role.eq(role), "family"].unique())
              for role in ["pronoun_anchor", "compound_anchor"]}
    rows = []
    for item in TARGETS:
        d = jaccard(x.loc[item].to_numpy(), aa.to_numpy())[0]
        for draw in range(2000):
            means = {}
            for role, families in groups.items():
                # Retain all forms of each sampled family, with multiplicity.
                indices = np.concatenate([np.flatnonzero(mm.family.eq(f).to_numpy())
                                          for f in rng.choice(families, len(families), replace=True)])
                means[role] = d[indices].mean()
            rows.append(dict(lemma=item, draw=draw+1, d_pronoun=means['pronoun_anchor'],
                             d_compound=means['compound_anchor'],
                             delta=means['pronoun_anchor']-means['compound_anchor']))
    pd.DataFrame(rows).to_csv(out / "comparator_resampling.csv", index=False)


def illustrate(x, manifest, blocks, excluded, out):
    features = ["compound_word", "paradigm_has_distinct_acc_form", "person_deixis", "reciprocal", "appears_in_subject"]
    examples = x.loc[["she", "herself", "everyone", *TARGETS], features]
    examples.to_csv(out / "illustrative_matrix.csv")
    dictionary = pd.DataFrame({"feature":x.columns, "block":blocks.values,
                               "presences":x.sum().values,
                               "interpretation":"1 = may exhibit the property; 0 = never does, in the source coding"})
    dictionary.to_csv(out / "feature_dictionary.csv", index=False)
    (out / "feature_filter.json").write_text(json.dumps({"excluded":excluded, "retained":len(x.columns),
        "block_counts":blocks.value_counts().to_dict(),
        "note":"Original column identifiers preserved, including spelling errors. Definitions are inherited coding claims, not independent grammatical validation."},indent=2)+"\n")


def plots(result, primary, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # Use the workspace house style when present; keep a standalone fallback.
    style = ROOT.parents[2] / ".house-style"
    if (style / "plot_style.py").exists():
        sys.path.insert(0, str(style))
        from plot_style import setup
        setup(font_size=11, tick_size=10)
    else:
        plt.rcParams.update({"font.family":"serif", "axes.spines.top":False,
                             "axes.spines.right":False, "legend.frameon":False})
    destination = out / "figures"
    destination.mkdir(exist_ok=True)
    def save(fig, name):
        for extension in ["pdf", "png"]:
            fig.savefig(destination / f"{name}.{extension}", dpi=300, bbox_inches="tight")
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 5.2))
    for role, label, color, marker in [
        ("pronoun_anchor", "Personal-pronoun controls", "#2E5077", "o"),
        ("compound_anchor", "Compound-determinative controls", "#C74F41", "s"),
        ("target", "Reciprocals", "#2D2D2D", "^")]:
        subset = primary.loc[primary.role.eq(role)]
        ax.scatter(subset.d_pronoun, subset.d_compound, label=label, c=color,
                   marker=marker, s=70 if role == "target" else 30, alpha=0.85)
    lo = min(primary.d_pronoun.min(), primary.d_compound.min()) - 0.035
    hi = max(primary.d_pronoun.max(), primary.d_compound.max()) + 0.045
    ax.plot([lo,hi],[lo,hi],color="0.55",lw=0.8,ls="--",label="Equal mean distances")
    for item, shift in [("she",(-10,30)),("everyone",(12,10)),("each_other",(-82,-18)),("one_another",(18,16))]:
        row = primary.loc[primary.lemma.eq(item)].iloc[0]
        ax.annotate(item.replace("_"," "),(row.d_pronoun,row.d_compound),xytext=shift,
                    textcoords="offset points",fontsize=10,
                    bbox=dict(facecolor="white",edgecolor="none",pad=1),
                    arrowprops=dict(arrowstyle="-",color="0.4",lw=.6))
    ax.set(xlim=(lo,hi),ylim=(lo,hi),xlabel="Mean Jaccard distance to personal-pronoun anchors",
           ylabel="Mean Jaccard distance to compound-determinative anchors")
    ax.set_aspect("equal",adjustable="box")
    ax.legend(loc="upper left", bbox_to_anchor=(0,1.25), fontsize=9)
    save(fig,"control_distances")
    base = result.query("pool == 'full' and holdout == 'item' and reference == 'fixed_total'")
    fig, axes = plt.subplots(1,2,figsize=(7.2,4.6),sharey=True)
    order = ["none","morph","synt","sem","phon"]
    colors = ["#2E5077","#C74F41","#3D825D","#946697"]
    for ax,item in zip(axes,TARGETS):
        for k,(metric,color) in enumerate(zip(METRICS,colors)):
            z = base.loc[base.lemma.eq(item)&base.metric.eq(metric)].set_index("drop").loc[order]
            label = {"jaccard":"Jaccard", "dice":"Dice", "hamming":"Hamming", "idf_jaccard":"Weighted Jaccard"}[metric]
            ax.scatter(np.arange(5)+(k-1.5)*0.12,z.delta_percentile,c=color,marker=["o","s","^","D"][k],label=label,s=35)
        ax.axhline(.5,c="0.6",ls="--",lw=.7)
        ax.set(xticks=np.arange(5),xticklabels=["All","No morph.","No syntax","No sem.","No phon."],
               ylim=(-.04,1.04), title=item.replace("_"," "))
        ax.tick_params(axis="x",rotation=30)
    axes[0].set_ylabel("Reference percentile of Δ")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=10)
    fig.tight_layout(rect=(0,.13,1,1))
    save(fig,"specification_percentiles")
    ref = pd.read_csv(out / "target_reference_draws.csv.gz")
    fig,axes = plt.subplots(2,3,figsize=(7.2,5.8),sharex="col",sharey=True)
    reference_fixed = ref.query("reference == 'fixed_total'")
    bins = {key: np.linspace(reference_fixed[key].min(), reference_fixed[key].max(), 36)
            for key in ["d_pronoun", "d_compound", "delta"]}
    for row_i,item in enumerate(TARGETS):
        draws = ref.query("lemma == @item and reference == 'fixed_total'")
        obs = primary.loc[primary.lemma.eq(item)].iloc[0]
        for ax,key,label in zip(axes[row_i],["d_pronoun","d_compound","delta"],
                               ["Distance to\npersonal pronouns","Distance to\ncompound determinatives","Δ (pronoun − compound)"]):
            ax.hist(draws[key],bins=bins[key],color="#E8E8E8",edgecolor="white")
            ax.axvline(obs[key],color="#2E5077",lw=1.5,label="Observed")
            ax.set_xlabel(label,fontsize=9)
        axes[row_i,0].set_ylabel(item.replace("_"," ")+"\nReference draws")
    axes[0,0].legend(fontsize=9)
    fig.tight_layout()
    save(fig,"target_references")


def control_summary(result, out):
    rows = []
    for key,z in result.query("reference == 'fixed_total'").groupby(["pool","holdout","drop","metric"]):
        c = z[z.eligible_control]
        correct = (c.role.eq("pronoun_anchor") & c.delta.lt(0)) | (c.role.eq("compound_anchor") & c.delta.gt(0))
        p_max = c.loc[c.role.eq("pronoun_anchor"),"delta"].max()
        c_min = c.loc[c.role.eq("compound_anchor"),"delta"].min()
        t = z[z.role.eq("target")]
        rows.append(dict(zip(["pool","holdout","drop","metric"],key)) |
                    dict(controls_correct=int(correct.sum()),n_controls=len(c),p_max=p_max,c_min=c_min,
                         both_targets_in_gap=bool((t.delta.gt(p_max)&t.delta.lt(c_min)).all())))
    pd.DataFrame(rows).to_csv(out/"control_recovery_by_specification.csv",index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws",type=int,default=None,help="Override only for smoke tests; recorded in metadata")
    parser.add_argument("--output",type=Path,default=ROOT / "data/revision/results")
    args = parser.parse_args()
    protocol,x,manifest,blocks,excluded = load_design()
    draws = args.draws if args.draws is not None else protocol["draws"]
    if draws < 2:
        parser.error("At least two reference draws required")
    out = args.output.resolve()
    out.mkdir(parents=True,exist_ok=True)
    illustrate(x,manifest,blocks,excluded,out)
    result,primary = analyse(x,manifest,blocks,protocol,draws,out)
    control_summary(result,out)
    coding_sensitivity(x,manifest,out)
    resample_families(x,manifest,protocol["seed"],out)
    plots(result,primary,out)
    meta = dict(protocol_sha256=hashlib.sha256((ROOT/"data/revision/protocol.json").read_bytes()).hexdigest(),
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                python=sys.version,numpy=np.__version__,pandas=pd.__version__,draws=draws,
                smoke_test=draws != protocol["draws"],rows=len(result),inputs=protocol["inputs"])
    (out/"run_metadata.json").write_text(json.dumps(meta,indent=2)+"\n")
    print(primary.loc[primary.lemma.isin([*TARGETS,"she","everyone"]),
          ["lemma","d_pronoun","d_compound","delta","d_pronoun_lower_tail","d_compound_lower_tail","delta_percentile"]].to_string(index=False))
    print(f"Wrote {len(result)} item/specification rows to {out}")


if __name__ == "__main__":
    main()
