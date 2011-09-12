#!/bwf/g64home/mbelhorn/.local/bin/python3
from array import array
import os
import fnmatch
import errno
import sys
import subprocess
import argparse

# Parse CLI for arguments and set default variables... FOR SCIENCE!
parser = argparse.ArgumentParser(
    description=('A BASF interface for the anomalous '
        'dilepton charge asymmetry analysis.\n'
        'The default mode (no arguments) runs a diagnostic test of the '
        'analysis. By specifying the appropriate flags, the analysis can be '
        'run on most data sets, real or monte carlo.'),
    epilog='Report bugs to M. Belhorn (matt.belhorn@gmail.com)')
parser.add_argument('-e', default = 0,
    type = int, dest = 'experiment_number', choices = [0, 43, 53, 67, 69, 71],
    help = 'process specific experiment number')
parser.add_argument('-s', default = -1,
    type = int, dest = 'stream_number',
    help = 'process specific monte carlo stream')
parser.add_argument('-f', default = 19299,
    type = int, dest = 'fs',
    help = 'process monte carlo generated with fs')
parser.add_argument('-c', default = False,
    action = 'store_true', dest = 'continuum',
    help = 'process continuum data')
parser.add_argument('-g', default = False,
    action = 'store_true', dest = 'generic_mc',
    help = 'process generic monte carlo')
parser.add_argument('-b', default = False,
    action = 'store_true', dest = 'jpsi_veto_os_only',
    help = 'allow charge bias in lepton pair vetoes')
parser.add_argument('-v', default = 0,
    action = 'count', dest = 'verbose_log',
    help = 'increase verbosity in log for each occurance')
parser.add_argument('-t', default = False,
    action = 'store_true', dest = 'test_run',
    help = 'test the code without running a BASF job')
options = parser.parse_args()

process_type = [
    options.experiment_number != 0,
    options.stream_number >= 0,
    options.continuum,
    options.generic_mc]

# Set the directories and instantiate processing variables.
analysis_directory = '/bwf/g64home/mbelhorn/analysis/adcab'
output_directory = analysis_directory + '/output'
mdst_directory = '/bwf/g64home/mbelhorn/analysis/MC/signalmc-u5s-mix/mdst'
process_files = []
output_name = 'Adcab.OOPS'

# Given inputs, find the data files to process.
if process_type == [0, 0, 0, 0]:
  print('Processing diagnostic run.')
  diagnostic_mdst = mdst_directory + '/fs19299.s0.e43.mBdBd.n0016545.mdst'
  options.stream_number = 0
  options.experiment_number = 43
  if os.path.isfile(diagnostic_mdst):
    print('on ' + diagnostic_mdst)
    process_files.append(diagnostic_mdst + ' 800')
  else:
    print('except mdst does not exist!')
  output_name = 'Adcab.diagnostic'

elif process_type == [1, 0, 0, 0]:
  print('Processing real data.')
  output_name = ('Adcab.DATA.e' + str(options.experiment_number))
  
elif process_type == [1, 1, 0, 0]:
  print('Processing signal MC.')
  mdst_base_filename = (
      'fs' + str(options.fs) + '.s' + str(options.stream_number) + 
      '.e' + str(options.experiment_number) + '.m*.n*.mdst')
  for file in os.listdir(mdst_directory):
    if fnmatch.fnmatch(file, mdst_base_filename):
      print('Found:', file)
      process_files.append(mdst_directory + '/' + str(file))
  output_name = ('Adcab.SGMC.fs' + str(options.fs)
      + '.s' + str(options.stream_number)
      + '.e' + str(options.experiment_number))
  if len(process_files) != 4:
    print('Warning: Number of MC files might not be right.')

elif process_type == [1, 0, 1, 0]:
  print('Processing continuum.')
  output_name = ('Adcab.SGMC.udsc.e' + str(options.experiment_number))

elif ((process_type == [1, 1, 1, 0]) or (process_type == [1, 1, 1, 1])):
  print('Processing generic MC continuum.')
  for event_type in ('charm', 'uds'):
    mdst_base_filename = ('dilepskim-e*' + str(options.experiment_number) + 
        'r*r*-s*' + str(options.stream_number) + '-evtgen-' +
        event_type + '-5S_onresonance-b*.index')
    mdst_directory = ('/bwf/g89home/g0mc/skim/skim5S/index/' + 
        'dilepskim/5S_onresonance/e0000' + str(options.experiment_number) + 
        '/evtgen-' + event_type + '/s0' + str(options.stream_number))
    for file in os.listdir(mdst_directory):
      if fnmatch.fnmatch(file, mdst_base_filename):
        print('Found:', file)
        process_files.append(mdst_directory + '/' + str(file))
  output_name = ('Adcab.GNMC.udsc.s' + str(options.stream_number) +
      '.e' + str(options.experiment_number))

elif process_type == [1, 1, 0, 1]:
  print('Processing generic MC B-Bbar events.')
  for event_type in ('bsbs', 'nonbsbs'):
    mdst_base_filename = ('dilepskim-e*' + str(options.experiment_number) + 
        'r*r*-s*' + str(options.stream_number) + '-evtgen-' +
        event_type + '-5S_onresonance-b*.index')
    mdst_directory = ('/bwf/g89home/g0mc/skim/skim5S/index/' + 
        'dilepskim/5S_onresonance/e0000' + str(options.experiment_number) + 
        '/evtgen-' + event_type + '/s0' + str(options.stream_number))
    for file in os.listdir(mdst_directory):
      if fnmatch.fnmatch(file, mdst_base_filename):
        print('Found:', file)
        process_files.append(mdst_directory + '/' + str(file))
  output_name = ('Adcab.GNMC.fs19299.s' + str(options.stream_number) +
      '.e' + str(options.experiment_number))

else:
  print('Input could not be interpereted.')
  print(options)
  sys.exit()

# If there are no files to be processed, we have a problem.
if len(process_files) == 0:
  print('No files to be processed. Check input.')
  sys.exit()

# Prepare the list of BASF initialization commands to use the Adcab module.
analysis_path_commands = [
    'path create main',
    'module register fix_mdst',
    'path add_module main fix_mdst',
    'path create analysis',
    'module register Adcab',
    'path add_module analysis Adcab',
    'module inquire_parameter Adcab',
    'path add_condition main <=:0:KILL',
    'path add_condition main >:0:analysis']

# Prepare the list of Adcab BASF parameters. BASF requires these to be strings
#   of the expected type, so we type-cast them to be sure.
analysis_parameters = [
    ['JPsi_Veto_OS_Only', str(int(options.jpsi_veto_os_only))],
    ['Verbose_Log',       str(int(options.verbose_log))],
    ['MC_Stream_Number',  str(int(options.stream_number))]]

# Test runs must exit before calling BASF.
if options.test_run:
  sys.exit()

# Set the log and BASF output paths.
# Warning! No safety check yet on whether these dirs exist!
basf_log = open(output_directory + '/logs/' + output_name + '.log', 'w') 
histogram_name = output_directory + '/hbks/' + output_name + '.hbk'

# Fire up BASF, with proper logs and a pipe to send it commands.
basf = subprocess.Popen('basf', stdin=subprocess.PIPE, stdout = basf_log,
    stderr = subprocess.STDOUT, shell = True)
basf_in = basf.stdin

# Send the initialization commands.
for command in analysis_path_commands:
  basf_in.write((command + '\n').encode('utf-8'))
  basf_in.flush()

# Send the module parameter settings.
for parameter in analysis_parameters:
  basf_in.write(
      ('module put_parameter Adcab ' + parameter[0] + '\\' + 
      parameter[1]+ '\n').encode('utf-8'))
  basf_in.flush()

# Finalize BASF's initialization.
basf_in.write(('initialize\n').encode('utf-8'))
basf_in.flush()
basf_in.write(('histogram define ' + histogram_name + '\n').encode('utf-8'))
basf_in.flush()

# Send the list of data files to process.
for file in process_files:
  basf_in.write(('process_event ' + file + '\n').encode('utf-8'))
  basf_in.flush()

# Close BASF and close the stream to the log. We're done here.
basf_in.write(('terminate\n').encode('utf-8'))
basf_in.flush()
basf_log.close()