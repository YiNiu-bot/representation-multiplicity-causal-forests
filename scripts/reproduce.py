#!/usr/bin/env python3
"""Refit the paper's experiments without overwriting included results."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def execute(command, env=None):
    print('Running:', ' '.join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), cwd=ROOT, env=env, check=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--stage', choices=['monte-carlo', 'empirical', 'all'], required=True)
    p.add_argument('--rscript', default='Rscript')
    p.add_argument('--empirical-rscript', default='Rscript')
    p.add_argument('--package-root', type=Path)
    p.add_argument('--replications', type=int, default=100)
    p.add_argument('--workers', type=int, default=3)
    p.add_argument('--output', type=Path, default=ROOT / 'rerun_outputs')
    a = p.parse_args()
    if a.replications < 1 or a.workers < 1:
        p.error('Replication and worker counts must be positive')
    output = a.output.resolve()
    if output == ROOT or (ROOT in output.parents and ROOT / 'rerun_outputs' not in [output, *output.parents]):
        p.error('Inside the repository, use a directory under rerun_outputs')
    output.mkdir(parents=True, exist_ok=True)
    if a.stage in ('monte-carlo', 'all'):
        code = ROOT / 'replication/code/monte_carlo/empirical'
        for name, script, prefix in [('dimension', 'run_default_mtry_dimension.R', 'RI_DIM'),
                                     ('seed_control', 'run_sign_seed_benchmark.R', 'RI_SEED')]:
            env = dict(os.environ)
            env.update({prefix + '_REPS': str(a.replications), 'RI_OUTPUT_DIR': str(output / name)})
            execute([a.rscript, code / script], env)
    if a.stage in ('empirical', 'all'):
        if a.package_root is None:
            p.error('--package-root is required for empirical refits')
        package = a.package_root.resolve()
        for name in ['primitives.R', 'stage0.R', 'cv.R', 'specifications_intersection.R']:
            if not (package / name).is_file():
                p.error('Missing official source file: ' + str(package / name))
        code = ROOT / 'replication/code/auto_dml'
        native, seed, canonical = [output / 'auto_dml' / name for name in ('native', 'seed_noise', 'canonical')]
        common = [f'--package_root={package}', '--tuning=theoretical', '--datasets=NSW,PSID,CPS',
                  '--specs=1,2', '--ntree=1000']
        execute([a.empirical_rscript, code / 'run_native.R', *common,
                 f'--output_dir={native}', '--phase=paired', '--seeds=1'])
        (seed / 'cache').mkdir(parents=True, exist_ok=True)
        for path in (native / 'cache').glob('*.rds'):
            shutil.copy2(path, seed / 'cache' / path.name)
        execute([a.empirical_rscript, code / 'run_native.R', *common,
                 f'--output_dir={seed}', '--phase=seed_noise', '--seeds=1,2,3,4,5'])
        execute([a.empirical_rscript, code / 'run_canonical_rank.R', native / 'cache', canonical, a.workers])


if __name__ == '__main__':
    main()
