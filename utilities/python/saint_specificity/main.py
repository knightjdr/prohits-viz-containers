import argparse
import math
import numpy as np
import pandas as pd

from distutils.util import strtobool

'''
Usage:

python3 main.py \
-c true \
-m fe \
-s saint.txt

output: saint-specificity.txt
'''

def specificity():
  options = parse_args()
  saint = read_saint(options.saint, options.control_subtract)
  calculate_specificity = get_specificty_calculator(saint, options.metric)
  saint_w_specificity = add_specificity_to_saint(saint, calculate_specificity)

  del saint_w_specificity['Abundance']
  del saint_w_specificity['Replicates']
  saint_w_specificity.to_csv('saint-specificity.txt', sep='\t', index=False)

def parse_args():
  parser = argparse.ArgumentParser(description='Perform GO enrichment')

  parser.add_argument(
    '--control_subtract', '-c',
    default=False,
    help='Subtract control average from spectral counts',
    type=lambda x: bool(strtobool(str(x)))
  )
  parser.add_argument(
    '--metric', '-m',
    default='fc',
    help='Specificity metric. Options: dscore, fe (default), sscore, wdscore and zscore',
  )
  parser.add_argument(
    '--saint', '-s',
    default='',
    help='SAINT file to process',
    required=True,
  )

  return parser.parse_args()

def read_saint(filename, control_subtract):
  '''
  Read is a SAINT file in tsv format and optionally subtract the average value
  in controls from the AvgSpec and the Spec column.
  '''
  df = pd.read_csv(filename, sep='\t')

  if control_subtract:
    df['Abundance'] = df.AvgSpec - df.ctrlCounts.str.split('|', expand=True).astype('float').mean(axis=1)
    df['Abundance'] = df['Abundance'].mask(df['Abundance'].lt(0), 0)
    df['Abundance'] = df['Abundance'].round(2)

    def subtract_control_from_reps(rep_str, ctrl_rep_str):
      reps = [float(x) for x in rep_str.split('|')]
      ctrl_mean = np.mean([float(x) for x in ctrl_rep_str.split('|')])
      return '|'.join(str(round(max(rep - ctrl_mean, 0), 2)) for rep in reps)
    df['Replicates'] = df.apply(lambda x: subtract_control_from_reps(x.Spec, x.ctrlCounts), axis = 1)
  else:
    df['Abundance'] = df.AvgSpec
    df['Replicates'] = df.Spec

  return df

def add_specificity_to_saint(saint, calculate_specificity):
  '''
  Add a specificity column to a SAINT file.
  '''
  df = saint
  df["Specificity"] = df.apply(lambda x: calculate_specificity(x.Prey, x.Abundance, x.Replicates), axis = 1)
  return df

def get_specificty_calculator(df, metric):
  '''
  Return a function for calculating the specificity for a given bait-prey pair.

  Options include: dscore, fe, sscore, wdscore and zscore.
  '''
  no_conditions = len(df.Bait.drop_duplicates())
  preys = df.Prey.drop_duplicates()
  spec_per_prey = {}
  for prey in preys:
    spec_per_prey[prey] = df.Abundance[df.Prey == prey].tolist()

  def dscore(prey, spec, reps):
    freq = no_conditions / len(spec_per_prey[prey])
    reproducibility = sum(float(x) > 0 for x in reps.split('|'))
    multiplier = pow(freq, reproducibility)
    adjusted_abundance = multiplier * spec
    return round(math.sqrt(adjusted_abundance), 2)

  def fe(prey, spec, *_):
    if spec == 0:
      return 0
    values_for_other_baits = spec_per_prey[prey].copy()
    values_for_other_baits.remove(spec)
    values_for_other_baits = np.pad(values_for_other_baits, (0, no_conditions - len(values_for_other_baits) - 1), 'constant')
    mean = np.mean(values_for_other_baits)
    return math.inf if mean == 0 else round(spec / mean, 2)

  def sscore(prey, spec, *_):
    freq = no_conditions / len(spec_per_prey[prey])
    adjusted_abundance = freq * spec
    return round(math.sqrt(adjusted_abundance), 2)
  
  def wdscore(prey, spec, reps):
    values_for_all_baits = spec_per_prey[prey].copy()
    values_for_all_baits = np.pad(values_for_all_baits, (0, no_conditions - len(values_for_all_baits)), 'constant')
    mean = np.mean(values_for_all_baits)
    sd = np.std(values_for_all_baits, ddof=1)
    omega = 1 if math.isnan(sd) or mean == 0 else sd / mean
    omega = max(omega, 1)
    freq = no_conditions / len(spec_per_prey[prey])
    weighted_frequency = freq * omega
    reproducibility = sum(float(x) > 0 for x in reps.split('|'))
    multiplier = pow(weighted_frequency, reproducibility)
    adjusted_abundance = multiplier * spec
    return round(math.sqrt(adjusted_abundance), 2)

  def zscore(prey, spec, *_):
    values_for_all_baits = spec_per_prey[prey].copy()
    values_for_all_baits = np.pad(values_for_all_baits, (0, no_conditions - len(values_for_all_baits)), 'constant')
    mean = np.mean(values_for_all_baits)
    sd = np.std(values_for_all_baits, ddof=1)
    return 0 if sd == 0 else round((spec - mean)/ sd, 2)

  if metric == 'dscore':
    return dscore
  if metric == 'sscore':
    return sscore
  if metric == 'wdscore':
    return wdscore
  if metric == 'zscore':
    return zscore
  return fe

if __name__ == "__main__":
  specificity()