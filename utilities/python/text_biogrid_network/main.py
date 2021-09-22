import argparse
import json
import pandas as pd
import re
import requests

from distutils.util import strtobool

'''
Usage:

python3 main.py \
-f file.txt \
-g genemap.json \
-k $access_key

output: biogrid-interactions.txt and cytoscape.txt
'''

def get_interactions():
  options = parse_args()

  ids, input_interactions = read_identifiers(options)
  id_map = read_id_map(options)
  mapped_ids, mapped_interactions = map_identifiers(ids, input_interactions, id_map)

  post_data = create_post_data(options, mapped_ids.keys())
  biogrid_data = fetch_interactions(post_data)
  write_biogrid_data(biogrid_data)

  interactions = extract_interaction_pairs(biogrid_data, mapped_ids.keys())
  interactions = merge_input_interactions(interactions, mapped_interactions)
  interactions = consolidate_symbols(options.id_type, interactions, mapped_ids, mapped_interactions)
  write_interactions(interactions, mapped_ids, options)

def parse_args():
  parser = argparse.ArgumentParser(description='Create a list of interactions retrieved from BioGRID from an input list')

  parser.add_argument(
    '--access_key', '-k',
    help='BioGRID access key',
    required=True,
  )
  parser.add_argument(
    '--evidence_list', '-el',
    default='',
    help='Pipe-separated list of evidence codes. Any interaction evidence with its '
      'Experimental System in the list will be excluded from the results unless '
      'include_evidence is set to true.',
  )
  parser.add_argument(
    '-fdr',
    default=0.01,
    help='FDR for filtering interactions from SAINT. Only used if both is_saint and '
    'include_saint_interactions are true',
    type=float
  )
  parser.add_argument(
    '--file', '-f',
    help='File with a newline, space or comma separated list of genes, or a SAINT file.',
    required=True,
  )
  parser.add_argument(
    '--id_type', '-it',
    default='symbol',
    help='Identifier type used in the input file.',
  )
  parser.add_argument(
    '--include_evidence', '-ie',
    default=True,
    help='If set to true, any interaction evidence with its Experimental System in '
      'the evidence_list will be included in the result',
    type=lambda x: bool(strtobool(str(x)))
  )
  parser.add_argument(
    '--include_primary_interactions', '-ipi',
    default=True,
    help='If true, in addition to interactions between genes on the gene_list, '
      ' interactions will also be fetched which have only one interactor on the '
      'gene_list i.e. the gene_list’s first order interactors will be included',
    type=lambda x: bool(strtobool(str(x)))
  )
  parser.add_argument(
    '--include_secondary_interactions', '-isi',
    default=False,
    help='If true interactions between the gene_list’s first order interactors '
      ' will be included. Ignored if include_interactors is false.',
    type=lambda x: bool(strtobool(str(x)))
  )
  parser.add_argument(
    '--include_saint_interactions', '-isai',
    default=True,
    help='If the input file is a SAINT file (is_saint is true), then include SAINT '
    'interactions in the output.',
    type=lambda x: bool(strtobool(str(x)))
  )
  parser.add_argument(
    '--inter_species_excluded', '-ise',
    default=False,
    help='If true, interactions with interactors from different species will be excluded.',
    type=lambda x: bool(strtobool(str(x)))
  )
  parser.add_argument(
    '--is_saint', '-is',
    default=False,
    help='Is the input file a SAINT formatted file.',
    type=lambda x: bool(strtobool(str(x)))
  )
  parser.add_argument(
    '--genemap', '-g',
    help='A file in JSON formatting mapping HUGO gene IDs to different identifiers',
    required=True,
  )
  parser.add_argument(
    '--max', '-m',
    default=10000,
    help='Number of results to fetch; this will be ignored if greater than 10,000, i.e. '
      'pagination using several requests is required to retrieve more than 10,000 interactions.',
    type=int
  )
  parser.add_argument(
    '--self_interactions_excluded', '-sie',
    default=True,
    help='If true, interactions with one interactor will be excluded.',
    type=lambda x: bool(strtobool(str(x)))
  )
  parser.add_argument(
    '--tax_id', '-ti',
    default='All',
    help='Pipe-separated list of NCBI taxonomy identifiers or All. Only genes from these '
    'organisms will be searched with reference to gene identifiers or names.',
  )
  parser.add_argument(
    '--throughput_tag', '-tt',
    default='any',
    help='If set to low or high, only interactions with Low throughput or High throughput '
      'in the throughput field will be returned. Interactions with both Low throughput and '
      'High throughput will be returned by either value.',
  )

  return parser.parse_args()

def read_identifiers(options):
  '''
  Read in a list of genes IDs from a file. The file can either be a list of genes
  separated by whitespace or comma, or a SAINT file. For a SAINT file the "Bait"
  column is read is as the genes, parsing the gene name as everything before the
  first underscore. Interactors from a SAINT file can also be returned. 
  '''

  file = options.file
  genes = []
  interactions = {}

  if options.is_saint:
    columns = ['Bait', 'PreyGene', 'BFDR']
    fdr = options.fdr
    include_saint_interactions = options.include_saint_interactions

    df = pd.read_csv(file, sep='\t', usecols=columns)
    df.Bait = df.Bait.str.split('_').str[0]
    genes = list(df.Bait.unique())

    if include_saint_interactions:
      df = df[df.BFDR <= fdr]
      interactions = {k: list(v) for k,v in df.groupby('Bait')['PreyGene']}
  else:
    with open(file, 'r') as f:
      text = f.read()
      genes = re.compile('[\s,]+').split(text.strip())
      genes = sorted(list(set(genes)))

  return genes, interactions

def read_id_map(options):
  '''
  Read in a file with a mapping of identifier types, and create a dict
  of the input identifier type to entrez. When creating a map from
  symbols, first add official symbols to the map, then aliases and finally
  previous symbols.
  '''

  genemapfile = options.genemap
  idtype = options.id_type

  with open(genemapfile) as json_data:
    data = json.load(json_data)

    genemap = {}
    string_ids = ['ensemblg', 'entrez']
    if idtype != 'symbol' and idtype not in string_ids:
      for ids in data.values():
        for id in ids[idtype]:
          genemap[id] = ids['entrez']

    if idtype == 'symbol':
      aliases = {}
      prev_symbols = {}
      symbols = {}
      for ids in data.values():
        symbols[ids['symbol']] = ids['entrez']
        for id in ids['aliasSymbol']:
          aliases[id] = ids['entrez']
        for id in ids['prevSymbol']:
          prev_symbols[id] = ids['entrez']
      genemap = {**prev_symbols, **aliases, **symbols}
  
    if idtype in string_ids:
      for ids in data.values():
        genemap[ids[idtype]] = ids['entrez']
    
    return genemap

'''
Case insensitive dictionary.
'''
class CaseInsensitiveDict(dict):
  @classmethod
  def _k(cls, key):
    return key.lower() if isinstance(key, str) else key

  def __init__(self, *args, **kwargs):
    super(CaseInsensitiveDict, self).__init__(*args, **kwargs)
    self._convert_keys()
  def __getitem__(self, key):
    return super(CaseInsensitiveDict, self).__getitem__(self.__class__._k(key))
  def __setitem__(self, key, value):
    super(CaseInsensitiveDict, self).__setitem__(self.__class__._k(key), value)
  def __delitem__(self, key):
    return super(CaseInsensitiveDict, self).__delitem__(self.__class__._k(key))
  def __contains__(self, key):
    return super(CaseInsensitiveDict, self).__contains__(self.__class__._k(key))
  def has_key(self, key):
    return super(CaseInsensitiveDict, self).has_key(self.__class__._k(key))
  def pop(self, key, *args, **kwargs):
    return super(CaseInsensitiveDict, self).pop(self.__class__._k(key), *args, **kwargs)
  def get(self, key, *args, **kwargs):
    return super(CaseInsensitiveDict, self).get(self.__class__._k(key), *args, **kwargs)
  def setdefault(self, key, *args, **kwargs):
    return super(CaseInsensitiveDict, self).setdefault(self.__class__._k(key), *args, **kwargs)
  def update(self, E={}, **F):
    super(CaseInsensitiveDict, self).update(self.__class__(E))
    super(CaseInsensitiveDict, self).update(self.__class__(**F))
  def _convert_keys(self):
    for k in list(self.keys()):
      v = super(CaseInsensitiveDict, self).pop(k)
      self.__setitem__(k, v)
    
def map_identifiers(ids, interactions, id_map):
  '''
  Create a mapping of input IDs to Entrez IDs. If prey interactions were extracted
  from a SAINT file, convert bait and prey IDs as well. Tests for identifier mapping
  are case insensitive.
  '''
  ci_id_map = CaseInsensitiveDict(id_map)
  mapped_ids = {ci_id_map[id.lower()]: id for id in ids if id.lower() in ci_id_map}

  mapped_interactions = {
      (ci_id_map[source.lower()], source): {ci_id_map[id.lower()]: { 'symbol': id } for id in targets if id.lower() in ci_id_map}
      for source, targets in interactions.items() if source.lower() in ci_id_map
  }
  return mapped_ids, mapped_interactions
  
def create_post_data(options, ids):
  '''
  Create data object for the POST request to BioGRID.
  '''

  data = {
    'accessKey': options.access_key,
    'additionalIdentifierTypes': 'ENTREZ_GENE',
    'format': 'json',
    'geneList': '|'.join(ids),
    'includeInteractors': options.include_primary_interactions,
    'includeInteractorInteractions': options.include_secondary_interactions,
    'interSpeciesExcluded': options.inter_species_excluded,
    'max': options.max,
    'selfInteractionsExcluded': options.self_interactions_excluded,
    'taxId': options.tax_id,
    'throughputTag': options.throughput_tag,
  }

  if options.evidence_list != '':
    data['evidenceList'] = options.evidence_list
    data['includeEvidence'] = options.include_evidence

  return data

def fetch_interactions(data):
  '''
  Fetch interactions from BioGRID
  '''
  url = 'https://webservice.thebiogrid.org/interactions/'
  r = requests.post(url, data=data)
  return [entry for entry in r.json().values()]

def write_biogrid_data(interactions):
  '''
  Write BioGRID data as a tsv file, with one interaction per row.
  '''
  with open('./biogrid-interactions.txt', 'w') as f:
    fields = [key for key in interactions[0].keys()]
    f.write('\t'.join(fields))
    f.write('\n')

    for interaction in interactions:
      values = [str(value) for value in interaction.values()]
      f.write('\t'.join(values))
      f.write('\n')

def extract_interaction_pairs(biogrid_data, ids):
  '''
  Create a minimal representation of the BioGRID data as a 2D dict
  indexed by input id tuple (source_id;source_symbol) and target id.
  '''
  interactions = {}

  for datum in biogrid_data:
    source_id = datum['ENTREZ_GENE_A']
    source_symbol = datum['OFFICIAL_SYMBOL_A']
    target_id = datum['ENTREZ_GENE_B']
    target_symbol = datum['OFFICIAL_SYMBOL_B']
    if target_id in ids:
      source_key = (target_id, target_symbol)
      target_id = source_id
      target_symbol = source_symbol
    else:
      source_key = (source_id, source_symbol)
    
    if source_key not in interactions:
      interactions[source_key] = {}
    interactions[source_key][target_id] = { 'symbol': target_symbol }

  return interactions

def merge_input_interactions(interactions, input_interactions):
  '''
  Merge any input interactions with data from BioGRID.
  '''
  sources = list(set([*interactions.keys(), *input_interactions.keys()]))

  merged = { source: {} for source in sources}
  for source in sources:
    targets = {
      **interactions.get(source, {}),
      **input_interactions.get(source, {}),
    }
    for target, target_data in targets.items():
      merged[source][target] = {
        'isprey': target in input_interactions.get(source, {}),
        'known': target in interactions.get(source, {}),
        'symbol': target_data['symbol'],
      }

  return merged

def consolidate_symbols(id_type, interactions, input_ids, input_interactions):
  '''
  Input gene symbols and those returned from BioGRID can disagree, so always use
  input symbols whenever possible.
  '''

  if id_type != 'symbol':
    return interactions

  input_ids_to_symbol = {}
  for targets in input_interactions.values():
    id_to_symbol = {id: v['symbol'] for id, v in targets.items()}
    input_ids_to_symbol = {
      **input_ids_to_symbol,
      **id_to_symbol,
    }
  input_ids_to_symbol = {
    **input_ids_to_symbol,
    **input_ids,
  }

  consolidated = {}
  for source, targets in interactions.items():
    source_key = (source[0], input_ids_to_symbol.get(source[0], source[1]))
    consolidated[source_key] = {
        target: {
          **target_data,
          'symbol': input_ids_to_symbol.get(target, target_data['symbol'])
        }
        for target, target_data in targets.items()
    }
  return consolidated

def write_interactions(interactions, input_ids, options):
  '''
  Write interactions as a tsv file.
  '''
  id_type = options.id_type
  include_saint_interactions = options.include_saint_interactions
  is_saint = options.is_saint

  with open('./cytoscape.txt', 'w') as f:
    if is_saint and include_saint_interactions:
      f.write('source\tsource Entrez\ttarget\ttarget Entrez\tis target a prey\tis target known\n')
      for source, targets in interactions.items():
        for target, target_data in targets.items():
          f.write(f'{source[1]}\t{source[0]}\t{target_data["symbol"]}\t{target}\t{target_data["isprey"]}\t{target_data["known"]}\n')
    elif id_type == 'symbol':
      f.write('source\tsource Entrez\ttarget\ttarget Entrez\n')
      for source, targets in interactions.items():
        for target, target_data in targets.items():
          f.write(f'{source[1]}\t{source[0]}\t{target_data["symbol"]}\t{target}\n')
    else:
      f.write('source ID\tsource symbol\tsource Entrez\ttarget symbol\ttarget Entrez\n')
      for source, targets in interactions.items():
        for target, target_data in targets.items():
          f.write(f'{input_ids.get(source[0], "")}\t{source[1]}\t{source[0]}\t{target_data["symbol"]}\t{target}\n')

if __name__ == '__main__':
  get_interactions()
