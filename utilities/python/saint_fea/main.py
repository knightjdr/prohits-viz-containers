import argparse
import csv
import os
import pandas as pd
from gprofiler import GProfiler
from pathlib import Path

'''
Usage:

python3 main.py \
-f 0.01 \
-s saint.txt \
-t 0

output: enrichment-saint.xlsx
'''

def enrich():
  options = parse_args()
  saint = read_saint(options)
  query, accessions_to_symbol = create_query_lists(saint)
  enrichment = gProfile(query, accessions_to_symbol)
  write_enrichment(enrichment, options.saint, options.top_preys)

def parse_args():
  parser = argparse.ArgumentParser(description='Perform GO enrichment')

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
  parser.add_argument(
    '--top_preys', '-t',
    default=0,
    help='Only use top preys for enrichment (default: %(default)d)',
    type=int,
  )

  return parser.parse_args()


def read_saint(options):
  fdr = options.fdr
  saintfile = options.saint
  top_preys = options.top_preys

  columns = ['Bait', 'Prey', 'PreyGene', 'AvgSpec', 'BFDR']
  df = pd.read_csv(saintfile, sep='\t', usecols=columns)
  df = df[df.BFDR <= fdr]

  if top_preys > 0:
    df = df.sort_values(['Bait', 'AvgSpec'], ascending=[True, False]).groupby(by='Bait').head(top_preys)

  df.reset_index(drop=True, inplace=True)

  return df

def create_query_lists(df):
  query = df.groupby('Bait')['Prey'].apply(list).to_dict()
  accessions_to_symbol = pd.Series(df.PreyGene.values, index=df.Prey).to_dict()
  return query, accessions_to_symbol

def gProfile(query, accessions_to_symbol):
  gp = GProfiler(return_dataframe=True)

  profile = gp.profile(
    domain_scope='annotated',
    no_evidences=False,
    no_iea=True,
    organism='hsapiens',
    query=query,
    significance_threshold_method='g_SCS',
    sources=['GO:MF','GO:CC','GO:BP','REAC','CORUM'],
    user_threshold=0.01,
  )

  def convert_accession_to_symbols(accesions):
    genes = [accessions_to_symbol[accession] for accession in accesions]
    genes.sort()
    return genes

  profile['genes'] = [convert_accession_to_symbols(x) for x in profile['intersections']]

  return profile

def write_enrichment(df, saintfile, top_preys):
  basename = os.path.basename(saintfile)
  filename = os.path.splitext(basename)[0]

  columns = [
    'bait',
    'source',
    'native',
    'name',
    'p_value',
    'term_size',
    'query_size',
    'intersection_size',
    'background_size',
    'precision',
    'recall',
    'genes',
  ]
  df = df.rename(columns={ 'query': 'bait', 'effective_domain_size': 'background_size' })
  df = df[columns]
  df = df.sort_values(['bait', 'p_value'], ascending=[True, True]).reset_index(drop=True)

  # Create excel file with all data
  BP = df[df.source == 'GO:BP']
  CC = df[df.source == 'GO:CC']
  MF = df[df.source == 'GO:MF']
  CORUM = df[df.source == 'CORUM']
  REAC = df[df.source == 'REAC']

  outfile = f'enrichment-{filename}.xlsx'
  if top_preys > 0:
    outfile = f'enrichment-top{top_preys}-{filename}.xlsx'

  # pylint: disable=abstract-class-instantiated
  with pd.ExcelWriter(outfile) as writer:
    df.to_excel(writer, index=False, sheet_name='all')
    BP.to_excel(writer, index=False, sheet_name='BP')
    CC.to_excel(writer, index=False, sheet_name='CC')
    MF.to_excel(writer, index=False, sheet_name='MF')
    CORUM.to_excel(writer, index=False, sheet_name='Corum')
    REAC.to_excel(writer, index=False, sheet_name='Reactome')

if __name__ == "__main__":
  enrich()
