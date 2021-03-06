import argparse

def parse_args():
  parser = argparse.ArgumentParser(description='Perform domain/motif enrichment')

  parser.add_argument(
    '--background', '-b',
    default='all',
    help='The background for the enrichment can be all proteins with annotated motifs (all [default]), '
      'or only proteins/genes found in the SAINT file (file)',
  )
  parser.add_argument(
    '--domains', '-d',
    default='',
    help='Domains for every gene in JSON format, with HUGO gene IDs as keys and an array of domain names',
    required=True,
  )
  parser.add_argument(
    '--fdr', '-f',
    default=0.01,
    help='FDR for significant preys (default: %(default).2f)',
    type=float,
  )
  parser.add_argument(
    '--genemap', '-g',
    default='',
    help='A file in JSON formatting mapping HUGO gene IDs to different identifiers',
    required=True,
  )
  parser.add_argument(
    '--idtype', '-i',
    default='refseqp',
    help='The type of identifier in the "Prey" column',
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
