import argparse

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
