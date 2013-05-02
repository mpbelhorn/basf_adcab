#!/usr/bin/python3
import subprocess

experiments = [31, 33, 35, 37, 41, 43, 45, 47, 49, 51, 53, 55, 61, 63, 65, 67, 69, 71, 73]
queue = 'b_a'

for experiment in experiments:
  lsf_args = [
      'bsub', '-q', queue,
      './run_basf.py',
      '-c',
      '-e', str(experiment)]
  log_args = lsf_args[3:]
  log_args.append('-t')
  log_process = subprocess.Popen(log_args)
  log_process.wait()
  lsf_process = subprocess.Popen(lsf_args)
  lsf_process.wait()
