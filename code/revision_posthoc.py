#!/usr/bin/env python3
"""Post hoc deterministic diagnostics for the revised reciprocal manuscript.

No imports from the original analysis pipeline and no random draws. Preserves all
original and recorded-follow-up outputs; writes only data/revision/posthoc/.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/revision/posthoc'
protocol_path = ROOT / 'data/revision/posthoc_protocol.json'
protocol = json.loads(protocol_path.read_text())
for name, digest in protocol['inputs'].items():
    if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
        raise ValueError('Input differs from post hoc reporting record: ' + name)
OUT.mkdir(exist_ok=True)
raw = pd.read_csv(ROOT / 'data/matrix_clean.csv').set_index('lemma').drop(columns='class')
x = raw.loc[:, raw.sum().ge(2)]
blocks = pd.read_csv(ROOT / 'data/feature_blocks_by_order.csv').set_index('feature').block.loc[x.columns]
manifest = pd.read_csv(ROOT / 'data/revision/comparators.csv').set_index('lemma')
active = manifest.index[manifest.role.ne('context')]
targets = ['each_other', 'one_another']
special = ['compound_word', 'does_not_inflect', 'has_derivationtionally_related_words']
components = protocol['component_columns']
assert set(components) <= set(blocks[blocks.eq('morph')].index)
assert components[-1] == 'yours' and len(components) == 49
other_morph = blocks[blocks.eq('morph')].index.difference(components).tolist()
remaining_morph = [c for c in other_morph if c not in special]
assert len(other_morph) == 16 and len(remaining_morph) == 13
assert other_morph == protocol['other16_morph_columns']
assert remaining_morph == protocol['remaining13_morph_columns']

def evaluate(name, cols, weighting='equal_feature', altered=False):
    frame = x.loc[:, cols].copy()
    if altered:
        for c in ['appears_in_subject', 'Functions_as_subject_of_interrogative_tag']:
            if c in frame: frame.loc['it_dum', c] = 1
    weights = np.ones(len(cols))
    if weighting == 'equal_block_total':
        bb = blocks.loc[cols]
        weights = np.array([1 / bb.eq(b).sum() for b in bb])
    rows = []
    for holdout in ['item', 'family']:
        for item in active:
            control = manifest.role.isin(['pronoun_anchor', 'compound_anchor'])
            if holdout == 'family' and manifest.loc[item, 'role'] != 'target':
                control &= manifest.family.ne(manifest.loc[item, 'family'])
            else:
                control.loc[item] = False
            a = frame.loc[control].to_numpy(dtype=bool)
            t = frame.loc[item].to_numpy(dtype=bool)
            numerator = np.sum(weights * np.logical_xor(a, t), axis=1)
            denominator = np.sum(weights * np.logical_or(a, t), axis=1)
            distances = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)
            pron = manifest.loc[control, 'role'].eq('pronoun_anchor').to_numpy()
            dp, dc = float(distances[pron].mean()), float(distances[~pron].mean())
            rows.append(dict(diagnostic=name, weighting=weighting, altered_dummy_it=altered,
                holdout=holdout, lemma=item, role=manifest.loc[item, 'role'],
                features=len(cols), d_pronoun=dp, d_compound=dc, delta=dp-dc))
    return pd.DataFrame(rows)

def without(*drop): return blocks.index[~blocks.isin(drop)].tolist()
specs = [
    ('all', x.columns.tolist(), 'equal_feature'),
    ('no_morph', without('morph'), 'equal_feature'),
    ('no_morph_no_phon', without('morph', 'phon'), 'equal_feature'),
    ('all', x.columns.tolist(), 'equal_block_total'),
    ('no_morph', without('morph'), 'equal_block_total'),
    ('no_morph_no_phon', without('morph', 'phon'), 'equal_block_total'),
    ('no_phon', without('phon'), 'equal_block_total'),
    ('no_morph_no_particle_stress', [c for c in without('morph') if c != 'Must_be_stressed_as_object_after_particle'], 'equal_block_total'),
    ('no_components', [c for c in x if c not in components], 'equal_feature'),
    ('no_other16_morph', [c for c in x if c not in other_morph], 'equal_feature'),
    ('no_shared3_morph', [c for c in x if c not in special], 'equal_feature'),
    ('no_remaining13_morph', [c for c in x if c not in remaining_morph], 'equal_feature'),
]
computed = pd.concat([evaluate(*spec) for spec in specs] +
    [evaluate(name, cols, weight, True) for name, cols, weight in specs if name in ['all', 'no_morph']], ignore_index=True)
computed.to_csv(OUT / 'deterministic_item_results.csv', index=False)
summary = []
for keys, g in computed.groupby(['diagnostic', 'weighting', 'altered_dummy_it', 'holdout']):
    p = g[g.role.eq('pronoun_anchor')].sort_values('delta', ascending=False)
    c = g[g.role.eq('compound_anchor')]
    bad = pd.concat([p[p.delta.ge(0)], c[c.delta.le(0)]])
    summary.append(dict(zip(['diagnostic', 'weighting', 'altered_dummy_it', 'holdout'], keys)) |
       dict(controls_recovered=57-len(bad), failures=bad.lemma.tolist(),
            pronoun_endpoint=p.iloc[0].lemma, pronoun_max=float(p.iloc[0].delta),
            next_pronoun=p.iloc[1].lemma, next_pronoun_delta=float(p.iloc[1].delta),
            compound_min=float(c.delta.min()),
            targets=g[g.role.eq('target')][['lemma', 'd_pronoun', 'd_compound', 'delta']].to_dict('records')))
saved = pd.read_csv(ROOT / 'data/revision/results/all_specifications.csv')
observed = saved[saved.reference.eq('fixed_total')]
endpoints = []
for keys, g in observed.groupby(['pool', 'holdout', 'drop', 'metric']):
    p = g[g.role.eq('pronoun_anchor') & g.eligible_control].sort_values('delta', ascending=False)
    c = g[g.role.eq('compound_anchor') & g.eligible_control]
    failures = pd.concat([p[p.delta.ge(0)], c[c.delta.le(0)]])
    endpoints.append(dict(zip(['pool', 'holdout', 'drop', 'metric'], keys)) |
      dict(pronoun_endpoint=p.iloc[0].lemma, pronoun_max=float(p.iloc[0].delta),
           next_pronoun=p.iloc[1].lemma, next_pronoun_delta=float(p.iloc[1].delta),
           failures=failures.lemma.tolist(),
           targets=g[g.role.eq('target')][['lemma', 'delta']].to_dict('records')))
assert len(endpoints) == 80
primary = observed.query("pool=='full' and holdout=='item' and drop=='none' and metric=='jaccard'")
benchmark = []
for role, other in [('pronoun_anchor', 'd_compound'), ('compound_anchor', 'd_pronoun')]:
    g = primary[primary.role.eq(role)]
    tail = g[other+'_lower_tail']
    benchmark.append(dict(role=role, non_designated_distance=other, n=len(g),
       median=float(tail.median()), maximum=float(tail.max()), below_005=int(tail.lt(.05).sum()),
       zero_hits=g.loc[np.isclose(tail, 1/20001, rtol=0, atol=1e-12), 'lemma'].tolist(),
       minimum_distance_to_other=float(g[other].min()),
       maximum_distance_to_own=float(g['d_pronoun' if role=='pronoun_anchor' else 'd_compound'].max()),
       maximum_distance_to_own_without_dummy=float(g[g.lemma.ne('it_dum')]['d_pronoun' if role=='pronoun_anchor' else 'd_compound'].max())))
feature_counts=[]
for col in special + blocks[blocks.eq('phon')].index.tolist() + ['inflects_for_case']:
    feature_counts.append(dict(feature=col, pronoun=int(x.loc[manifest.role.eq('pronoun_anchor'), col].sum()),
        compound=int(x.loc[manifest.role.eq('compound_anchor'), col].sum()),
        targets=x.loc[targets, col].tolist()))
# Verify the independent implementation against the existing saved baseline rows.
errors=[]
for name, drop in [('all', 'none'), ('no_morph', 'morph')]:
    for weight in ['equal_feature', 'equal_block_total']:
        a=computed.query('diagnostic==@name and weighting==@weight and not altered_dummy_it')
        b=(observed.query("pool=='full' and metric=='jaccard' and drop==@drop") if weight=='equal_feature'
           else pd.read_csv(ROOT/'data/revision/followup/block_weight_item_results.csv').query('drop==@drop'))
        z=a.merge(b,on=['lemma','holdout'],suffixes=('_check','_saved'))
        error=max(float((z[c+'_check']-z[c+'_saved']).abs().max()) for c in ['d_pronoun','d_compound','delta'])
        assert len(z)==118 and error<1e-12
        errors.append(dict(diagnostic=name,weighting=weight,rows=len(z),maximum_absolute_error=error))
result=dict(kind='post_hoc_deterministic_review_checks', new_random_draws=0,
    component_columns=components, other16_morph_columns=other_morph, remaining13_morph_columns=remaining_morph,
    baseline_crosschecks=errors, diagnostics=summary, original_grid_endpoints=endpoints,
    control_benchmark=benchmark, feature_counts=feature_counts,
    inputs={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in [
      'data/matrix_clean.csv','data/feature_blocks_by_order.csv','data/revision/comparators.csv',
      'data/revision/results/all_specifications.csv','data/revision/followup/block_weight_item_results.csv']},
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(OUT/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
summary_frame = pd.json_normalize(summary)
summary_frame.to_csv(OUT/'diagnostic_summary.csv',index=False)
pd.DataFrame(feature_counts).to_csv(OUT/'feature_counts.csv',index=False)
pd.DataFrame(benchmark).to_csv(OUT/'control_benchmark.csv',index=False)
pd.DataFrame(endpoints).to_csv(OUT/'original_grid_endpoints.csv',index=False)
metadata = dict(protocol_sha256=hashlib.sha256(protocol_path.read_bytes()).hexdigest(),
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    inputs=protocol['inputs'], numpy=np.__version__, pandas=pd.__version__,
    item_result_rows=len(computed), diagnostic_summary_rows=len(summary),
    original_observed_specifications=len(endpoints), new_random_draws=0,
    baseline_crosschecks=errors, outputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
        for p in OUT.iterdir() if p.is_file() and p.name!='metadata.json'})
(OUT/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
for r in summary:
    if r['holdout']=='item':
        print(r['diagnostic'],r['weighting'],'altered='+str(r['altered_dummy_it']),
              [round(t['delta'],8) for t in r['targets']],r['controls_recovered'],r['failures'])
print('ALL ORIGINAL ENDPOINTS:',sorted(set(r['pronoun_endpoint'] for r in endpoints)))
print('ALL ORIGINAL FAILURES:',sorted(set(i for r in endpoints for i in r['failures'])))
print('CONTROL BENCHMARK:',json.dumps(benchmark))
print('FEATURE COUNTS:',json.dumps(feature_counts))
print('REMAINING 13:',remaining_morph)
print('BASELINE ERRORS:',json.dumps(errors))
