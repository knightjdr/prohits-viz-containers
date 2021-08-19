import argparse
import os
import pandas as pd

'''
Usage:

python3 main.py \
-f folder \
-t bagel

output: [tool]-converted.txt
'''

def convert():
  options = parse_args()

  files = get_files(options.folder)
  
  df = pd.DataFrame()
  if options.tool == 'bagel':
    df = merge_and_add_conditions(files)
  if options.tool == 'drugz':
    df = merge_and_add_conditions(files)
  if options.tool == 'mageck':
    df = convert_mageck(files)
  if options.tool == 'ranks':
    df = convert_ranks(files)

  df.to_csv(f'{options.tool}-converted.txt', sep='\t', index=False)

def parse_args():
  parser = argparse.ArgumentParser(description='Convert CRISPR output files to a single file compatible with ProHits-viz')

  parser.add_argument(
    '--folder', '-f',
    help='Folder containing the files to merge/convert',
    required=True,
  )
  parser.add_argument(
    '--tool', '-t',
    help='The tool used for CRISPR analysis. Should be one of "bagel", "drugz", "mageck" or "ranks"',
    required=True,
  )

  return parser.parse_args()

def get_files(folder):
  files = os.listdir(folder)
  return [f'{folder}/{file}' for file in files]

def extract_condition_from_filename(file):
  base = os.path.basename(file)
  return os.path.splitext(base)[0]

def get_column_names(file):
  df = pd.read_csv(file, sep='\t')
  return df.columns.values.tolist()

def merge_and_add_conditions(files):
  '''
  Files are simply merged together, adding one addition column to the start
  specifying the condition (using the filename)
  '''
  data = []
  for file in files:
    df = pd.read_csv(file, sep='\t')

    condition = extract_condition_from_filename(file)
    df['condition'] = condition
    data.append(df)

  merged = pd.concat(data, axis=0)
  merged.reset_index(drop=True, inplace=True)
  return (move_condition_column(merged))

def move_condition_column(df):
  return df[ ['condition'] + [ col for col in df.columns if col != 'condition' ] ]

def convert_mageck(files):
  '''
  MAGeCK has two output formats. The first from the "test" command, or RRA, just requires
  merging of files and adding the condition name. The second from the "mle" command has each
  condition with its own set of columns in a single file. For example:

  | Gene  | sgRNA | HL60|beta | HL60|z | HL60|p-value | HL60|fdr | HL60|wald-p-value | HL60|wald-fdr | KBM7|beta | etc... |
  '''
  columns = get_column_names(files[0])

  # Check for "test" output format
  if 'neg|score' in columns:
    return merge_and_add_conditions(files)
  
  def extract_condition(column):
    return column.split('|')[0]
  
  def filter_condition_columns(column):
    return 'beta' in column

  data = []
  desired_column_names = ['Gene', 'sgRNA', 'beta', 'z', 'p-value', 'fdr', 'wald-p-value', 'wald-fdr']
  for file in files:
    df = pd.read_csv(file, sep='\t')
    columns = get_column_names(file)
    conditions = list(map(extract_condition, filter(filter_condition_columns, columns)))

    for condition in conditions:
      condition_columns = ['Gene', 'sgRNA', f'{condition}|beta', f'{condition}|z', f'{condition}|p-value', f'{condition}|fdr', f'{condition}|wald-p-value', f'{condition}|wald-fdr']
      df_partial = df[condition_columns].copy()
      df_partial.columns = desired_column_names
      df_partial['condition'] = condition
      data.append(df_partial)

  merged = pd.concat(data, axis=0)
  merged.reset_index(drop=True, inplace=True)
  return (move_condition_column(merged))

def convert_ranks(files):
  '''
  RANKS has a single output format, but for v1 it can have two entries for the same gene, one
  with a negative (depletion) score and one with a positive score. Sometimes it may have only one score if
  only depletion was calculated. v2 produces a single score for each gene encapsulating both
  depletion and enrichment.
  '''
  merged = merge_and_add_conditions(files)

  # Check for v1 formatted data will duplicated rows
  if merged.duplicated(['condition','Gene']).any():
    df_depletion = merged[merged['RANKS_score'] <= 0]
    df_depletion.columns = ['condition', 'Gene', 'depletion_score', 'depletion_p-value', 'depletion_FDR', 'depletion_#_of_sgRNAs_considered']
    df_enrichment = merged[merged['RANKS_score'] > 0]
    df_enrichment.columns = ['condition', 'Gene', 'enrichment_score', 'enrichment_p-value', 'enrichment_FDR', 'enrichment_#_of_sgRNAs_considered']
    df = pd.merge(df_depletion, df_enrichment, how='outer', on=['condition', 'Gene'])
  
    nanFillValues = {
      'depletion_score': 0,
      'depletion_p-value': 1,
      'depletion_FDR': 1,
      'depletion_#_of_sgRNAs_considered': 0,
      'enrichment_score': 0,
      'enrichment_p-value': 1,
      'enrichment_FDR': 1,
      'enrichment_#_of_sgRNAs_considered': 0,
    }
    df.fillna(value=nanFillValues, inplace=True)
    df[['depletion_#_of_sgRNAs_considered', 'enrichment_#_of_sgRNAs_considered']] = df[['depletion_#_of_sgRNAs_considered', 'enrichment_#_of_sgRNAs_considered']].astype('int64')
    return df;

  return merged

if __name__ == "__main__":
  convert()