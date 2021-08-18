import argparse
import csv
import statistics

'''
Usage:

python3 main.py \
-f 0.01 \
-s saint.txt

output: saint-statistics.txt
'''

def summarize():
  options = parse_args()
  summary = read_saint(options.saint, options.fdr)
  write_summary(summary, options.saint)

def parse_args():
  parser = argparse.ArgumentParser(description='Calculate interaction stats for SAINT file')

  parser.add_argument(
    '--fdr', '-f',
    default=0.01,
    help='FDR for significant preys (default: %(default).2f)',
    type=float,
  )
  parser.add_argument(
    '--saint', '-s',
    default='',
    help='SAINT file to process',
    required=True,
  )

  return parser.parse_args()

def read_saint(filename, cutoff):
  summary = {
    'interactions': {
      'bait': {},
      'total': {},
      'unique': {},
    },
    'preys': {
      'significant': {},
      'total': {},
    },
  }

  with open(filename) as csv_file:
    reader = csv.reader(csv_file, delimiter = '\t')
    header = next(reader, None)
    fdr_column_index = header.index('BFDR')

    for row in reader:
      bait = row[0]
      bait_gene = bait.split('_')[0]
      prey_gene = row[2]
      fdr = float(row[fdr_column_index])

      
      summary['preys']['total'][prey_gene] = True

      if fdr <= cutoff:
        genes = [bait_gene, prey_gene]
        genes.sort()
        summary['interactions']['total'][f'{bait}-{prey_gene}'] = True
        summary['interactions']['unique'][f'{genes[0]}-{genes[1]}'] = True
        summary['preys']['significant'][prey_gene] = True

        if bait not in summary['interactions']['bait']:
          summary['interactions']['bait'][bait] = 0
        summary['interactions']['bait'][bait] += 1


  return summary

def write_summary(summary, saint_file):
  with open('saint-statistics.txt', 'w') as summary_file:
    summary_file.write(f'file: {saint_file}\n')
    summary_file.write(f'interactions: unique - {len(summary["interactions"]["unique"])}, total - {len(summary["interactions"]["total"])}\n')
    summary_file.write(f'preys: significant - {len(summary["preys"]["significant"])}, total - {len(summary["preys"]["total"])}\n')
    
    min_bait, min_preys = get_min_prey_number(summary)
    summary_file.write(f'min: {min_bait} with {min_preys} prey(s)\n')

    max_bait, max_preys = get_max_prey_number(summary)
    summary_file.write(f'max: {max_bait} with {max_preys} prey(s)\n')

    summary_file.write(f'median: {get_median_prey_number(summary)} prey(s)\n')

def get_max_prey_number(summary):
  max_preys = -1
  max_bait = ''
  for bait, total in summary['interactions']['bait'].items():
    if total > max_preys:
      max_bait = bait
      max_preys = total
    
  return max_bait, max_preys

def get_min_prey_number(summary):
  min_preys = -1
  min_bait = ''
  for bait, total in summary['interactions']['bait'].items():
    if min_preys == -1 or total < min_preys:
      min_bait = bait
      min_preys = total
    
  return min_bait, min_preys

def get_median_prey_number(summary):
  prey_counts = [total for _, total in summary['interactions']['bait'].items()]
  return statistics.median(prey_counts)

if __name__ == "__main__":
  summarize()
