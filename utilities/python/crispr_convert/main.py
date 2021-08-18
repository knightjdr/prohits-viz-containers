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
  
  if options.tool == 'bagel':
    df = convert_bagel(files)

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

def move_condition_column(df):
  return df[ ['condition'] + [ col for col in df.columns if col != 'condition' ] ]

def convert_bagel(files):
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

if __name__ == "__main__":
  convert()