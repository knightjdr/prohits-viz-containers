import argparse
import csv
import datetime
import os.path
import pandas as pd
import re

from openpyxl import load_workbook
from pathlib import Path

'''
Usage:

python3 main.py \
-f file.txt \
-c column1|column2

output: file-fixed.txt and summary.txt
'''

def fix_symbols():
  options = parse_args()
  
  is_excel = check_file_in_excel_format(options.file)
  data = read_file(options.file, is_excel)
  symbols_to_fix = read_symbols_to_fix()

  fixed, summary = fix_data(data, options.columns, symbols_to_fix)
  write_summary(summary, is_excel, symbols_to_fix)
  write_fixed_file(fixed, options.file, is_excel)

def parse_args():
  parser = argparse.ArgumentParser(description='Fix gene symbols that have been converted to dates by excel')

  parser.add_argument(
    '--columns', '-c',
    default='',
    help='Pipe-separated list of columns to fix. It will fix all columns if none are supplied.',
  )
  parser.add_argument(
    '--file', '-f',
    help='Text or Excel file',
    required=True,
  )

  return parser.parse_args()

def check_file_in_excel_format(filepath):
  '''
  Checks if a file is in Excel format. Only .xlsx is acceptable.
  '''
  extension = os.path.splitext(filepath)[1]
  try:
    load_workbook(filepath)
  except Exception as inst:
    return False
  else:
    return True

def read_file(filepath, is_excel):
  '''
  Read a tab-delimited text or Excel file into a dict keyed
  by sheet name.
  '''
  data = {}

  if is_excel:
    data = pd.read_excel(filepath, sheet_name=None)
  else:
    data['sheet'] = pd.read_csv(filepath, sep='\t')

  return data

def read_symbols_to_fix():
  '''
  Read in a list of symbols that represent erroneous conversions by Excel.
  '''
  symbols_to_fix = {}

  p = Path(__file__).with_name('symbols-to-fix.txt')
  with p.open('r') as csv_file:
    reader = csv.reader(csv_file, delimiter='\t')
    header = next(reader, None)

    for row in reader:
      excel_symbol = row[0]
      converted_symbol = row[1]
      official_symbol = row[2]
      ambiguous = row[3] == 'True'
      if excel_symbol not in symbols_to_fix:
        symbols_to_fix[excel_symbol] = {
          'ambiguous': False,
          'converted_symbols': [],
          'official_symbols': [],
        }
      symbols_to_fix[excel_symbol]['ambiguous'] = ambiguous
      symbols_to_fix[excel_symbol]['converted_symbols'].append(converted_symbol)
      symbols_to_fix[excel_symbol]['official_symbols'].append(official_symbol)

  for key, data in symbols_to_fix.items():
    symbols_to_fix[key]['converted_symbols'] = sorted(list(set(symbols_to_fix[key]['converted_symbols'])))
    symbols_to_fix[key]['official_symbols'] = sorted(list(set(symbols_to_fix[key]['official_symbols'])))

  return symbols_to_fix

def fix_data(data, columns, symbols_to_fix):
  '''
  Check columns for symbols to fix and fix them whenever they can be unambiguously assigned
  to a single gene. As a first step in the process, anything that looks like
  a numerical date in the format of 2021-09-21 will get converted to its day-month representation
  (e.g. 21-Sep). Then symbols are tested. It will return the fixed dataframe and summary of the changes.
  '''
  columns_to_check = columns.split('|') if columns != '' else []
  unambiguous_map = {symbol: symbol_data['official_symbols'][0] for symbol, symbol_data in symbols_to_fix.items() if not symbol_data['ambiguous']}

  def replace_numerical_date(matchobj):
    (month, day) = matchobj.groups()
    month_abbr = datetime.date(1900, int(month), 1).strftime('%B')[0:3]
    return f'{day}-{month_abbr}'

  fixed = {}
  pattern_numerical_dates = re.compile('^20\d{2}-0*(\d{1,2})-0*(\d{1,2}).*')
  pattern_symbols_to_fix = '|'.join(symbols_to_fix.keys())
  summary = []
  for sheet, df in data.items():
    fixed[sheet] = df
    range = columns_to_check if len(columns_to_check) > 0 else df.columns.values.tolist()
    for column in range:
      df[column] = df[column].astype(str).apply(lambda a: pattern_numerical_dates.sub(replace_numerical_date, a))
      matches = df[df[column].str.contains(pattern_symbols_to_fix)][column]
      unique_matches = pd.unique(matches)
      column_summary = [(sheet, column, match) for match in pd.unique(matches)]
      summary.extend(column_summary)

      fixed[sheet][column] = df[column].map(unambiguous_map).fillna(df[column])

  return fixed, summary

def write_summary(summary, is_excel, symbols_to_fix):
  '''
  Write a summary of changed symbols.
  '''
  start, end = (0, 5) if is_excel else (1, 5)

  column_headings = ['sheet', 'column', 'original symbol', 'fixed', 'correct(ed) symbol(s)']

  with open('./summary.txt', 'w') as f:
    f.write('\t'.join(column_headings[start:end]))
    f.write('\n')

    for entry in summary:
      f.write('\t'.join(entry[start:3]))
      f.write(f'\t{not symbols_to_fix[entry[2]]["ambiguous"]}')
      f.write(f'\t{", ".join(symbols_to_fix[entry[2]]["official_symbols"])}')
      f.write('\n')

def write_fixed_file(data, filepath, is_excel):
  '''
  Write the original file with fixed symbols in the correct input format.
  '''
  base = os.path.basename(filepath)
  filename = os.path.splitext(base)[0]

  if is_excel:
    outfile = f'{filename}-fixed.xlsx'
    with pd.ExcelWriter(outfile) as writer:
      for sheet, df in data.items():
        df.to_excel(writer, index=False, sheet_name=sheet)
  else:
    outfile = f'{filename}-fixed.txt'
    data['sheet'].to_csv(outfile, sep='\t', index=False)

if __name__ == '__main__':
  fix_symbols()