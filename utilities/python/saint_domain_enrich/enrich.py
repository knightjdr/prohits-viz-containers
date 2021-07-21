import csv
import json
import os
import pandas as pd
import scipy.stats as stats

def enrich(options):
  domains = read_domains(options.domains)
  genemap = read_gene_map(options)
  
  saint = read_saint(options.saint)
  saint_mapped = map_file_ids(saint, genemap)
  filtered_saint = filter_saint(options, saint_mapped)

  background = get_background(options, saint_mapped, domains)

  domains_by_id, ids_by_domain = parse_domains(domains, background)
  domains_by_bait = count_domains_by_bait(filtered_saint, domains_by_id)
  enriched_domains_by_bait = calculate_enrichment(domains_by_bait, ids_by_domain, len(background), options.fdr)

  write_enriched_domains(options, enriched_domains_by_bait)

def read_domains(domainfile):
  with open(domainfile) as json_data:
    return json.load(json_data)

def read_gene_map(options):
  genemapfile = options.genemap
  idtype = options.idtype

  def parse_list_ids(ids, genemap, value):
    for id in ids[idtype]:
      genemap[id] = value
  def parse_string_id(ids, genemap, value):
    genemap[ids[idtype]] = value
  
  parse_ids = parse_list_ids
  if idtype == 'entrez':
    parse_ids = parse_string_id

  with open(genemapfile) as json_data:
    data = json.load(json_data)

    genemap = {}
    for hugoid, ids in data.items():
      parse_ids(ids, genemap, hugoid)
    return genemap

def read_saint(saintfile):
  columns = ['Bait', 'Prey', 'PreyGene', 'AvgSpec', 'BFDR']
  return pd.read_csv(saintfile, sep='\t', usecols=columns)

def map_file_ids(saint, genemap):
  mapped = saint
  mapped.Prey = mapped.Prey.str.split('.').str[0]
  mapped.Prey = mapped.Prey.map(genemap).fillna(mapped.Prey)
  return mapped

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

def parse_domains(domain_list_by_id, background_ids):
  ids_by_domain = {}
  domains_by_id = {}

  for id, domains in domain_list_by_id.items():
    if id in background_ids:
      elements = {}
      for domain in domains:
        if domain['name'] not in ids_by_domain:
          ids_by_domain[domain['name']] = []

        if domain['name'] not in elements:
          ids_by_domain[domain['name']].append(id)
          elements[domain['name']] = {
            'count': 0,
            'length': 0,
          }
        elements[domain['name']]['count'] += 1
        elements[domain['name']]['length'] += (domain['end'] - domain['start'] + 1)

      domains_by_id[id] = elements

  return domains_by_id, ids_by_domain

def count_domains_by_bait(saint, domains_by_id):
  domains_by_bait = {}

  baits = list(saint.Bait.unique())
  preyLookup = pd.Series(saint.PreyGene.values, index=saint.Prey).to_dict()

  for bait in baits:
    domains_by_bait[bait] = {
      'preysInDatabase': 0,
      'domain': {},
    }

    prey_ids = list(saint.Prey[saint.Bait == bait])
    for prey_id in prey_ids:
      if prey_id and prey_id in domains_by_id:
        domains_by_bait[bait]['preysInDatabase'] += 1

        for domain, domain_data in domains_by_id[prey_id].items():
          if domain not in domains_by_bait[bait]['domain']:
            domains_by_bait[bait]['domain'][domain] = {
              'countByAccession': [],
              'lengthByAccession': [],
              'preys': [],
            }

          domains_by_bait[bait]['domain'][domain]['countByAccession'].append(domain_data['count'])
          domains_by_bait[bait]['domain'][domain]['lengthByAccession'].append(domain_data['length'])
          domains_by_bait[bait]['domain'][domain]['preys'].append(preyLookup[prey_id])

  return domains_by_bait

def calculate_enrichment(domains_by_bait, ids_by_domain, background_size, fdr):
  enriched_domains_by_bait = {}

  def perform_domain_enrichment(bait_data, preys_in_database):
    # sourcery skip: inline-immediately-returned-variable, list-comprehension
    pvalues = {}
    stats = {}

    for domain, domain_data in bait_data.items():
      prey_count = len(domain_data['preys'])
      bait_fold_enrichment = prey_count / preys_in_database
      background_size_w_domain = len(ids_by_domain[domain])
      background_fold_enrichment = background_size_w_domain / background_size
      fold_enrichment = bait_fold_enrichment / background_fold_enrichment

      pvalue = fishers_test(
        prey_count,
        preys_in_database,
        background_size_w_domain,
        background_size,
      )
      pvalues[domain] = pvalue
      stats[domain] = {
        'background_size_w_domain': background_size_w_domain,
        'fold_enrichment': fold_enrichment,
        'pvalue': pvalue,
      }

    adj_pvalues, corrected_fdr = bh_correction(pvalues, fdr)

    enriched_domains = []
    for domain, adj_pvalue in sort_dict_by_value(adj_pvalues).items():
      if adj_pvalue == 0 or adj_pvalue < corrected_fdr[domain]:
        enriched_domains.append({
          'domain': domain,
          'no_genes_with_domain': len(bait_data[domain]['preys']),
          'no_genes': preys_in_database,
          'background_size_w_domain': stats[domain]['background_size_w_domain'],
          'background_size': background_size,
          'fold_enrichment': stats[domain]['fold_enrichment'],
          'pvalue': stats[domain]['pvalue'],
          'adj_pvalue': adj_pvalue,
          'bh_fdr': corrected_fdr[domain],
          'genes': bait_data[domain]['preys'],
        })
    
    return enriched_domains

  for bait, bait_data in domains_by_bait.items():
    enriched_domains_by_bait[bait] = perform_domain_enrichment(bait_data['domain'], bait_data['preysInDatabase'])

  return enriched_domains_by_bait

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
    adj_pvalue = min(adj_pvalue, 1)
    if last_adjusted < adj_pvalue:
      adj_pvalue = last_adjusted
    last_adjusted = adj_pvalue
    adjusted_pvalues[key] = adj_pvalue
    corrected_fdr[key] = fdr * ranks[key] / no_tests

  return adjusted_pvalues, corrected_fdr

def write_enriched_domains(options, enriched_domains_by_bait):
  saintfile = options.saint
  top_preys = options.top_preys

  basename = os.path.basename(saintfile)
  filename = os.path.splitext(basename)[0]

  outfile = f'domain-enrichment-{filename}.xlsx'
  if top_preys > 0:
    outfile = f'domain-enrichment-top{top_preys}-{filename}.xlsx'

  domains = []

  for bait, bait_data in enriched_domains_by_bait.items():
    for domain in bait_data:
      domains.append({
        'bait': bait,
        **domain,
      })

  # pylint: disable=abstract-class-instantiated
  with pd.ExcelWriter(outfile) as writer:
    pd.DataFrame(domains).to_excel(writer, index=False, sheet_name='domains')