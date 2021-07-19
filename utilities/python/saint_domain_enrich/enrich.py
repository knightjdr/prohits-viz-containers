import csv
import json
import os
import pandas as pd
import scipy.stats as stats

def enrich(options):
  domains = read_domains(options.domains)
  saint = read_saint(options.saint)
  filtered_saint = read_saint(options, saint)
  background = get_background(options, saint, domains)

  # sequence_elements_by_accession, sequence_elements = read_sequence_elements(options.gix, background)
  # elements_by_bait = count_elements_by_bait(saint, sequence_elements_by_accession)
  # enriched_elements_by_bait = calculate_enrichment(elements_by_bait, sequence_elements, len(background), options.fdr_enrichment)
  # write_enriched_elements(enriched_elements_by_bait, options.saint, options.top_preys)

def read_domains(domainfile):
  with open(domainfile) as json_data:
    return json.load(json_data)

def read_saint(saintfile):
  columns = ['Bait', 'Prey', 'PreyGene', 'AvgSpec', 'BFDR']
  return pd.read_csv(saintfile, sep='\t', usecols=columns)

def filter_saint(options, saint):
  fdr = options.fdr
  top_preys = options.top_preys

  df = saint[saint.BFDR <= fdr]

  if top_preys > 0:
    df = df.sort_values(['Bait', 'AvgSpec'], ascending=[True, False]).groupby(by='Bait').head(top_preys)

  df.reset_index(drop=True, inplace=True)

  return df

def get_background(options, saint, domains):
  background = options.background
  if background == 'all':
    return list(domains.keys())

  return list(saint.Prey.unique())

def read_sequence_elements(db_file, background):
  sequence_elements = {}
  sequence_elements_by_accession = {}

  with open(db_file) as json_data:
    data = json.load(json_data)

    for entry in data:
      if any(accession in background for accession in entry['uniprot']):
        accession = entry['uniprot'][0]
        symbol = entry['gene']
    
        elements = {}
        for element in entry['domains']:
          if element['name'] not in sequence_elements:
            sequence_elements[element['name']] = {}
          sequence_elements[element['name']][accession] = symbol

          if element['name'] not in elements:
            elements[element['name']] = {
              'count': 0,
              'length': 0,
              'type': element['type'],
            }
          elements[element['name']]['count'] += 1
          elements[element['name']]['length'] += (element['end'] - element['start'] + 1)

        for uniprotID in entry['uniprot']:
          sequence_elements_by_accession[uniprotID] = elements

  return sequence_elements_by_accession, sequence_elements

def count_elements_by_bait(saint, sequence_elements):
  elements_by_bait = {}

  baits = list(saint.Bait.unique())
  preyLookup = pd.Series(saint.PreyGene.values, index=saint.Prey).to_dict()

  for bait in baits:
    elements_by_bait[bait] = {
      'preysInDatabase': 0,
      'domain': {},
      'region': {},
    }

    accessions = list(saint.Prey[saint.Bait == bait])
    for accession in accessions:
      if accession and accession in sequence_elements:
        elements_by_bait[bait]['preysInDatabase'] += 1

        for element, element_data in sequence_elements[accession].items():
          element_type = element_data['type']

          if element not in elements_by_bait[bait][element_type]:
            elements_by_bait[bait][element_type][element] = {
              'countByAccession': [],
              'lengthByAccession': [],
              'preys': [],
            }

          elements_by_bait[bait][element_type][element]['countByAccession'].append(element_data['count'])
          elements_by_bait[bait][element_type][element]['lengthByAccession'].append(element_data['length'])
          elements_by_bait[bait][element_type][element]['preys'].append(preyLookup[accession])

  return elements_by_bait

def calculate_enrichment(elements_by_bait, sequence_elements, background_size, fdr):
  enriched_elements_by_bait = {}

  def perform_element_type_enrichment(bait_data, element_type):
    pvalues = {}
    stats = {}

    for element, element_data in bait_data[element_type].items():
      prey_count = len(element_data['preys'])
      bait_fold_enrichment = prey_count / bait_data['preysInDatabase']
      background_size_w_element = len(sequence_elements[element])
      background_fold_enrichment = background_size_w_element / background_size
      fold_enrichment = bait_fold_enrichment / background_fold_enrichment
  
      pvalue = fishers_test(
        prey_count,
        bait_data['preysInDatabase'],
        background_size_w_element,
        background_size,
      )
      pvalues[element] = pvalue
      stats[element] = {
        'background_size_w_element': background_size_w_element,
        'fold_enrichment': fold_enrichment,
        'pvalue': pvalue,
      }

    adj_pvalues, corrected_fdr = bh_correction(pvalues, fdr)

    enriched_elements = []
    for element, adj_pvalue in sort_dict_by_value(adj_pvalues).items():
      if adj_pvalue == 0 or adj_pvalue < corrected_fdr[element]:
        enriched_elements.append({
          'element': element,
          'no_genes_with_element': len(bait_data[element_type][element]['preys']),
          'no_genes': bait_data['preysInDatabase'],
          'background_size_w_element': stats[element]['background_size_w_element'],
          'background_size': background_size,
          'fold_enrichment': stats[element]['fold_enrichment'],
          'pvalue': stats[element]['pvalue'],
          'adj_pvalue': adj_pvalue,
          'bh_fdr': corrected_fdr[element],
          'genes': bait_data[element_type][element]['preys'],
        })
    return enriched_elements

  for bait, bait_data in elements_by_bait.items():
    enriched_elements_by_bait[bait] = {
      'domain': perform_element_type_enrichment(bait_data, 'domain'),
      'region': perform_element_type_enrichment(bait_data, 'region'),
    }

  return enriched_elements_by_bait

def sort_dict_by_value(dict):
  return {k: v for k, v in sorted(dict.items(), key=lambda item: item[1])}

def fishers_test(n11, n1p, np1, npp):
  '''
                bait  background_wo_bait  | background
  with_term     n11   n12                 | n1p
  without_term  n21   n22                 | n2p
                ---------------------------
                np1   np2                   npp
  '''
  n12 = n1p - n11
  n21 = np1 - n11
  n22 = npp - n1p - n21
  _, pvalue = stats.fisher_exact([[n11, n12], [n21, n22]], 'greater')
  return pvalue

def bh_correction(pvalues, fdr):
	# order p-values
  sorted_pvalues = sort_dict_by_value(pvalues)

  # determine rank/order for each key
  last_pvalue = -1
  rank = 0
  ranks = {}
  for key, pvalue in sorted_pvalues.items():
    if pvalue > last_pvalue :
      rank += 1
    ranks[key] = rank
    last_pvalue = pvalue

  # calculate adjusted p-values and fdr for each key
  no_tests = len(sorted_pvalues)
  adjusted_pvalues = {}
  corrected_fdr = {}
  last_adjusted = 1
  for key in reversed(list(sorted_pvalues.keys())):
    pvalue = sorted_pvalues[key]

    adj_pvalue = pvalue * no_tests / ranks[key]
    if adj_pvalue > 1:
      adj_pvalue = 1
    if last_adjusted < adj_pvalue:
      adj_pvalue = last_adjusted
    last_adjusted = adj_pvalue
    adjusted_pvalues[key] = adj_pvalue
    corrected_fdr[key] = fdr * ranks[key] / no_tests

  return adjusted_pvalues, corrected_fdr

def write_enriched_elements(enriched_elements, saintfile, top_preys):
  basename = os.path.basename(saintfile)
  filename = os.path.splitext(basename)[0]

  outfile = f'data/processed/dr-enrichment-{filename}.xlsx'
  if top_preys > 0:
    outfile = f'data/processed/dr-enrichment-top{top_preys}-{filename}.xlsx'

  domains = []
  regions = []

  for bait, bait_data in enriched_elements.items():
    for domain in bait_data['domain']:
      domains.append({
        'bait': bait,
        **domain,
      })
    for region in bait_data['region']:
      regions.append({
        'bait': bait,
        **region,
      })

  # pylint: disable=abstract-class-instantiated
  with pd.ExcelWriter(outfile) as writer:
    pd.DataFrame(domains).to_excel(writer, index=False, sheet_name='domains')
    pd.DataFrame(regions).to_excel(writer, index=False, sheet_name='regions')