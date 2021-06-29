import argparse

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
