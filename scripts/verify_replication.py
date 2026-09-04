#!/usr/bin/env python3
"""Verify frozen outputs and reconstruct the manuscript's numerical tables."""
import csv
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path
import statistics

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'replication/results'


def read(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def close(a, b, tol=1e-10):
    assert math.isfinite(float(a)) and math.isfinite(float(b))
    assert abs(float(a) - float(b)) <= tol, (a, b)


def summary(rows, keys):
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[k] for k in keys)].append(row)
    output = {}
    for key, values in groups.items():
        assert sorted(int(v['replication']) for v in values) == list(range(1, 101)), key
        x = [float(v['value']) for v in values]
        assert all(math.isfinite(v) for v in x)
        output[key] = {'mean': statistics.mean(x), 'mcse': statistics.stdev(x) / 10}
    return output


def main():
    manifest = read(ROOT / 'replication/metadata/manifest.csv')
    listed = set()
    for row in manifest:
        file = ROOT / row['file']
        assert file.is_file() and not file.is_symlink(), row['file']
        assert hashlib.sha256(file.read_bytes()).hexdigest() == row['sha256'], row['file']
        listed.add(row['file'])
    actual = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file()
              and '.git' not in p.parts and 'rerun_outputs' not in p.parts
              and '__pycache__' not in p.parts
              and p.suffix not in ('.aux', '.bbl', '.blg', '.fdb_latexmk', '.fls', '.log', '.out')
              and str(p.relative_to(ROOT)) not in ('replication/metadata/manifest.csv', 'paper/source/main.pdf')}
    assert actual == listed, {'unlisted': sorted(actual-listed), 'missing': sorted(listed-actual)}

    mc = RESULTS / 'monte_carlo'
    dimension = summary(read(mc / 'dimension_metrics.csv'), ['dimension', 'method', 'metric'])
    seed = summary(read(mc / 'seed_benchmark_metrics.csv'), ['comparison', 'metric'])
    losses_raw = read(mc / 'prediction_losses.csv')
    assert len(losses_raw) == 8000
    losses = summary(losses_raw, ['dimension', 'method', 'metric'])
    indexed = {(r['replication'], r['dimension'], r['method'], r['metric']): float(r['value']) for r in losses_raw}
    assert len(indexed) == 8000
    for r in range(1, 101):
        for metric in ('mse', 'mae', 'mean_error', 'wrong_sign', 'regret', 'value_gain', 'oracle_gain', 'treatment_share'):
            close(indexed[str(r), '500', 'class_x1', metric], indexed[str(r), '500', 'class_x2', metric], 0)
        for p in ('40', '100', '250', '500'):
            methods = ['ordinary_x1', 'ordinary_x2'] + (['class_x1', 'class_x2'] if p == '500' else [])
            for method in methods:
                close(indexed[str(r), p, method, 'regret'],
                      indexed[str(r), p, method, 'oracle_gain'] - indexed[str(r), p, method, 'value_gain'])
    pairs_raw = read(mc / 'paired_loss_differences.csv')
    for row in pairs_raw:
        a, b = row['method'].split('_minus_')
        r, p, m = row['replication'], row['dimension'], row['metric']
        close(row['value'], indexed[r, p, a, m] - indexed[r, p, b, m])
    paired = summary(pairs_raw, ['dimension', 'method', 'metric'])
    for name, calculated in [('mc_loss_summary.csv', losses), ('mc_paired_loss_summary.csv', paired)]:
        for row in read(mc / name):
            z = calculated[row['dimension'], row['method'], row['metric']]
            close(row['mean'], z['mean'])
            close(row['mcse'], z['mcse'])

    source = (ROOT / 'paper/source/main.tex').read_text()
    def displayed(number, places, commas=False):
        token = format(number, (',' if commas else '') + f'.{places}f')
        assert token in source, ('Number absent from manuscript source', token)

    def table_row(tokens):
        row = ' & '.join(tokens) + r'\\'
        assert row in source, ('Empirical table row does not match outputs', row)

    def estimate(row):
        return f"{float(row['atet']):,.1f} ({float(row['se']):,.1f})"

    table1 = []
    for p in ('40', '100', '250', '500'):
        row = {'dimension': int(p)}
        for metric in ('sign_reversal', 'strong_reversal_010'):
            z = dimension[p, 'ordinary_default', metric]
            row[metric] = {k: 100*v for k, v in z.items()}
            displayed(100*z['mean'], 2)
        table1.append(row)
    close(dimension['500', 'class_sampled_default', 'mean_absolute_prediction_change']['mean'], 0, 0)
    table1.append({'seed_control': seed['second_seed', 'sign_reversal']})
    for key, z in losses.items():
        if key[2] in ('mse', 'regret') and key[1] != 'class_x2':
            displayed(z['mean'], 4)
            displayed(z['mcse'], 4)
    for key, z in paired.items():
        if key[1].startswith('class_') and key[2] in ('mse', 'regret'):
            displayed(z['mean'], 4)
            displayed(z['mcse'], 4)

    empirical = RESULTS / 'empirical'
    native = read(empirical / 'native/paired_seed_results.csv')
    canonical = read(empirical / 'canonical/results.csv')
    assert len(native) == 18 and len(canonical) == 12
    seed_native = read(empirical / 'native_seed/paired_seed_results.csv')
    assert len(seed_native) == 30
    table3, table4 = [], []
    for dataset in ('NSW', 'PSID', 'CPS'):
        for spec in ('1', '2'):
            def choose(rows):
                return [r for r in rows if r['dataset'] == dataset and r['specification'] == spec]
            n = {r['variant']: r for r in choose(native)}
            c = {r['variant']: r for r in choose(canonical)}
            full = n['published_full']
            for row in list(n.values()) + list(c.values()):
                assert int(row['seed']) == 1 and int(row['ntree']) == 1000
                assert float(row['se']) > 0 and int(row['treated']) == 185
                close(row['ci_low'], float(row['atet'])-1.96*float(row['se']), 1e-7)
                close(row['ci_high'], float(row['atet'])+1.96*float(row['se']), 1e-7)
                displayed(float(row['atet']), 1, True)
                displayed(float(row['se']), 1, True)
            old_seed = {int(r['seed']): float(r['atet']) for r in choose(seed_native)}
            assert set(old_seed) == set(range(1, 6))
            close(old_seed[1], full['atet'], 1e-7)
            changes = [abs(old_seed[s] - old_seed[1]) for s in range(2, 6)]
            label = f'{dataset} / {spec}'
            table_row([label, estimate(full), estimate(n['quotient_default']), estimate(n['quotient_fixed_mtry'])])
            table_row([label,
                       f"{float(n['quotient_default']['atet'])-float(full['atet']):+.1f}",
                       f"{float(n['quotient_fixed_mtry']['atet'])-float(full['atet']):+.1f}",
                       f'{statistics.mean(changes):.1f} ({max(changes):.1f})'])
            table3.append({'dataset': dataset, 'specification': int(spec),
                           'native': n, 'seed_mean_absolute_change': statistics.mean(changes),
                           'seed_maximum_change': max(changes)})
            assert int(c['canonical_full']['mtry']) == int(c['canonical_reduced']['mtry']) == (8 if spec == '1' else 10)
            delta = float(c['canonical_reduced']['atet']) - float(c['canonical_full']['atet'])
            diag = next(r for r in choose(read(empirical / 'canonical/diagnostics.csv')) if r['subset'] == 'all')
            close(delta, diag['atet_difference'], 1e-7)
            table_row([label, estimate(c['canonical_full']), estimate(c['canonical_reduced']),
                       f'{delta:+.1f}', f"{100*float(diag['contrast_sign_disagreement']):.2f}"])
            table4.append({'dataset': dataset, 'specification': int(spec), 'canonical': c, 'change': delta})
            for row in read(empirical / 'canonical' / f'{dataset.lower()}_spec{spec}_certificate.csv'):
                assert row['identical_canonical_values'] == 'TRUE'
            for row in read(empirical / 'canonical' / f'{dataset.lower()}_spec{spec}_score_checks.csv'):
                close(row['score_reconstruction_gap'], 0, 1e-8)
    assert 'Chi, Chien-Ming' in (ROOT / 'paper/source/references.bib').read_text()
    print(json.dumps({'status': 'PASS', 'files_checked': len(manifest),
                      'scope': 'Frozen-output reconstruction, not forest refitting or proof validation',
                      'table1': table1,
                      'table2_losses': {' | '.join(k): v for k, v in losses.items() if k[2] in ('mse', 'regret')},
                      'table2_paired': {' | '.join(k): v for k, v in paired.items() if k[1].startswith('class_') and k[2] in ('mse', 'regret')},
                      'table3': table3, 'table4': table4}, indent=2))


if __name__ == '__main__':
    main()
