#!/bwf/g64home/mbelhorn/.local/bin/python3
import subprocess

experiments = [43, 53, 67, 69, 71]
streams = range(9)
queue = 'g'

for experiment in experiments:
  for stream in streams:
    lsf_args = [
        'bsub', '-q', queue,
        './run_basf.py',
        '-s', str(stream),
        '-e', str(experiment)]
    log_args = lsf_args[3:]
    log_args.append('-t')
    log_process = subprocess.Popen(log_args)
    log_process.wait()
    lsf_process = subprocess.Popen(lsf_args)
    lsf_process.wait()
