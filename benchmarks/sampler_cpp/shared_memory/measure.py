"""Fresh processes for timings; separate runs for peak RAM across all workers."""
import csv
import os
import platform
from pathlib import Path
import subprocess
import time

root = Path(os.getenv('EMC2_BENCHMARK_ROOT', '/home/nsteven/2026/StanCompare/emc2_shared'))
results = root / 'results'
results.mkdir(exist_ok=True)
macos = platform.system() == 'Darwin'
ram_column = 'peak_rss_mib' if macos else 'peak_pss_mib'


def pss_group(group):
    total = 0
    for path in Path('/proc').iterdir():
        if not path.name.isdigit():
            continue
        try:
            if os.getpgid(int(path.name)) == group:
                for line in (path / 'smaps_rollup').read_text().splitlines():
                    if line.startswith('Pss:'):
                        total += int(line.split()[1])
        except (ProcessLookupError, FileNotFoundError, PermissionError):
            pass
    return total


def rss_group(group):
    rows = subprocess.check_output(['ps', '-axo', 'pgid=,rss='], text=True)
    return sum(int(rss) for pgid, rss in (line.split() for line in rows.splitlines())
               if int(pgid) == group)


rows = []
for cores in ([2] if macos else [2, 15]):
    for repeat in range(4):
        memory = repeat == 3
        prefix = results / f'shared-{cores}-{repeat}'
        peak = 0
        with prefix.with_suffix('.log').open('w') as log:
            job = subprocess.Popen(['Rscript', '--vanilla', str(Path(__file__).with_name('benchmark.R')),
                                    str(cores), str(prefix)],
                                   stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            if memory:
                while job.poll() is None:
                    peak = max(peak, (rss_group if macos else pss_group)(job.pid))
                    time.sleep(.05 if macos else .02)
            else:
                job.wait()
        if job.returncode:
            raise RuntimeError(prefix.with_suffix('.log').read_text())
        with prefix.with_suffix('.csv').open() as f:
            row = next(csv.DictReader(f))
        row.update(repeat=repeat, memory_run=memory)
        row[ram_column] = peak / 1024 if memory else ''
        rows.append(row)
        print(row, flush=True)
        with (results / 'measurements.csv').open('w') as f:
            writer = csv.DictWriter(f, fieldnames=list(row))
            writer.writeheader()
            writer.writerows(rows)
