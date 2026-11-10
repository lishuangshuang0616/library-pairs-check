"""Compare two release binaries on generated run_macos.py data (macOS)."""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=Path, required=True)
parser.add_argument('--after', type=Path, required=True)
parser.add_argument('--data-root', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
work = Path(tempfile.mkdtemp(prefix='pair-optimization-'))
# Real sample directories and n*n comparisons; identical content is a load test,
# not biological validation. Hardlinks avoid copying the generated inputs.
for count in [4, 16]:
    for sample in range(count):
        for role, kind in [('cDNA', 'cdna'), ('oligo', 'oligo')]:
            dest = work / f'samples_{count}' / f'S{sample:03}' / role
            dest.mkdir(parents=True)
            os.link(args.data_root / f'small_plain_{kind}.fq', dest / 'lane_R1.fq')
scenarios = [(name, args.data_root / name, threads) for name in
             ['medium_gzip', 'large_gzip', 'singletons_gzip', 'four_by_four_gzip']
             for threads in [1, 8]]
scenarios += [(f'samples_{n}', work / f'samples_{n}', 8) for n in [4, 16]]
rows = []
for name, source, threads in scenarios:
    runs = {'before': [], 'after': []}
    for repeat in range(3):
        reports = {}
        for version in (['before', 'after'] if repeat % 2 == 0 else ['after', 'before']):
            output = work / f'{name}-{threads}-{version}.json'
            start = time.perf_counter()
            result = subprocess.run(['/usr/bin/time', '-l', str(getattr(args, version).resolve()),
                                     '--input', str(source), '--threads', str(threads),
                                     '--output', str(output)], capture_output=True, text=True, check=True)
            elapsed = time.perf_counter() - start
            report = json.loads(output.read_text())
            reports[version] = report
            runs[version].append({'wall_seconds': elapsed,
                                  'rss_mib': int(re.search(r'(\d+)\s+maximum resident set size', result.stderr)[1]) / 2**20,
                                  'compare_seconds': report['timings']['compare_seconds'],
                                  'rank_seconds': report['timings']['rank_seconds']})
        # Optimizations must preserve every numerical pair metric exactly.
        after_pairs = [{k: pair[k] for k in before_pair}
                       for before_pair, pair in zip(reports['before']['pairs'], reports['after']['pairs'])]
        assert reports['before']['pairs'] == after_pairs, name
    row = {'scenario': name, 'requested_threads': threads,
           'pairs': len(reports['after']['pairs']), 'runs': runs}
    for version in runs:
        row[version] = {'median_seconds': statistics.median(r['wall_seconds'] for r in runs[version]),
                        'max_rss_mib': max(r['rss_mib'] for r in runs[version])}
    rows.append(row)
    print(name, threads, row['pairs'], row['before'], row['after'], flush=True)
args.output.write_text(json.dumps({'work': str(work), 'data_root': str(args.data_root),
                                   'notes': 'Synthetic warm-cache data; 3 interleaved repeats; all numerical pair metrics identical.',
                                   'results': rows}, indent=2))
