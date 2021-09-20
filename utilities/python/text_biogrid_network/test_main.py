from os import access
import pyfakefs.fake_filesystem_unittest
import unittest

from .main import (
  consolidate_symbols,
  create_post_data,
  extract_interaction_pairs,
  map_identifiers,
  merge_input_interactions,
  read_id_map,
  read_identifiers,
  write_biogrid_data,
  write_interactions,
)

class ReadIdentifiers(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()

  def compare_actual_expected(self, options, expected_ids, expected_interactions):
    actual_ids, actual_interactions = read_identifiers(options)
    self.assertEqual(actual_ids, expected_ids)
    self.assertEqual(actual_interactions, expected_interactions)

  def test_list(self):
    file_contents = (
      'AAA\n'
      'BBB CCC   DDD\n'
      'EEE,FFF  GGG\n'
      'AAA,BBB\n'
    )
    filepath = '/test/file.txt'
    self.fs.create_file(filepath, contents=file_contents)

    class Options:
      file = filepath
      is_saint = False
    options = Options()

    expected_ids = ['AAA', 'BBB', 'CCC', 'DDD', 'EEE', 'FFF', 'GGG']
    expected_interactions = {}

    self.compare_actual_expected(options, expected_ids, expected_interactions)

  def test_saint_file(self):
    file_contents = (
      'Bait\tPrey\tPreyGene\tSpec\tAvgSpec\tBFDR\n'
      'AAA\tP11111\tprey1\t\t10\t0.01\n'
      'AAA\tP22222\tprey2\t\t20\t0\n'
      'AAA\tP33333\tprey3\t\t30\t0.02\n'
      'AAA\tP44444\tprey4\t\t15\t0.01\n'
      'AAA\tP55555\tprey5\t\t25\t0.01\n'
      'AAA\tP66666\tprey6\t\t40\t0.01\n'
      'BBB_Nt\tP11111\tprey1\t\t10\t0.05\n'
      'BBB_Nt\tP22222\tprey2\t\t20\t0.01\n'
      'BBB_Nt\tP33333\tprey3\t\t30\t0.01\n'
      'CCC\tP11111\tprey1\t\t10\t0.05\n'
    )
    filepath = '/test/file.txt'
    self.fs.create_file(filepath, contents=file_contents)

    class Options:
      fdr = 0.01
      file = filepath
      include_saint_interactions = False
      is_saint = True
    options = Options()

    expected_ids = ['AAA', 'BBB', 'CCC']
    expected_interactions = {}

    self.compare_actual_expected(options, expected_ids, expected_interactions)

  def test_saint_file_with_interactions(self):
    file_contents = (
      'Bait\tPrey\tPreyGene\tSpec\tAvgSpec\tBFDR\n'
      'AAA\tP11111\tprey1\t\t10\t0.01\n'
      'AAA\tP22222\tprey2\t\t20\t0\n'
      'AAA\tP33333\tprey3\t\t30\t0.02\n'
      'AAA\tP44444\tprey4\t\t15\t0.01\n'
      'AAA\tP55555\tprey5\t\t25\t0.01\n'
      'AAA\tP66666\tprey6\t\t40\t0.01\n'
      'BBB_Nt\tP11111\tprey1\t\t10\t0.05\n'
      'BBB_Nt\tP22222\tprey2\t\t20\t0.01\n'
      'BBB_Nt\tP33333\tprey3\t\t30\t0.01\n'
      'CCC\tP11111\tprey1\t\t10\t0.05\n'
    )
    filepath = '/test/file.txt'
    self.fs.create_file(filepath, contents=file_contents)

    class Options:
      fdr = 0.01
      file = filepath
      include_saint_interactions = True
      is_saint = True
    options = Options()

    expected_ids = ['AAA', 'BBB', 'CCC']
    expected_interactions = {
      'AAA': ['prey1', 'prey2', 'prey4', 'prey5', 'prey6'],
      'BBB': ['prey2', 'prey3'],
    }

    self.compare_actual_expected(options, expected_ids, expected_interactions)

class ReadIdMap(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()

  def get_test_options(self, arg_idtype):
    file_contents = (
      '{\n'
      '"1": {\n'
      '"aliasSymbol": ["aaa"],\n'
      '"ensemblg": "ENSG00000000111",\n'
      '"entrez": "11",\n'
      '"refseqp": ["NP_11111"],\n'
      '"prevSymbol": ["aa" ,"a1"],\n'
      '"symbol": "AAA"\n'
      '},\n'
      '"2": {\n'
      '"aliasSymbol": ["AAA"],\n'
      '"ensemblg": "ENSG00000000222",\n'
      '"entrez": "22",\n'
      '"refseqp": ["NP_22222", "NP_02222"],\n'
      '"prevSymbol": ["bb"],\n'
      '"symbol": "BBB"\n'
      '},\n'
      '"3": {\n'
      '"aliasSymbol": ["aa"],\n'
      '"ensemblg": "ENSG00000000333",\n'
      '"entrez": "33",\n'
      '"refseqp": ["NP_33333"],\n'
      '"prevSymbol": [],\n'
      '"symbol": "CCC"\n'
      '},\n'
      '"5": {\n'
      '"aliasSymbol": [],\n'
      '"ensemblg": "ENSG00000000555",\n'
      '"entrez": "55",\n'
      '"refseqp": ["NP_55555"],\n'
      '"prevSymbol": ["ee"],\n'
      '"symbol": "EEE"\n'
      '},\n'
      '"6": {\n'
      '"aliasSymbol": [],\n'
      '"ensemblg": "ENSG00000000666",\n'
      '"entrez": "66",\n'
      '"refseqp": ["NP_66666"],\n'
      '"prevSymbol": [],\n'
      '"symbol": "FFF"\n'
      '}\n'
      '}\n'
    )
    filepath = '/test/genemap.json'
    self.fs.create_file(filepath, contents=file_contents)

    class Options:
      genemap = filepath
      id_type = arg_idtype

    return Options()

  def test_list_ids(self):
    options = self.get_test_options('refseqp')
    expected = {
      'NP_11111': '11',
      'NP_22222': '22',
      'NP_02222': '22',
      'NP_33333': '33',
      'NP_55555': '55',
      'NP_66666': '66',
    }

    self.assertEqual(read_id_map(options), expected)

  def test_string_id(self):
    options = self.get_test_options('ensemblg')
    expected = {
      'ENSG00000000111': '11',
      'ENSG00000000222': '22',
      'ENSG00000000333': '33',
      'ENSG00000000555': '55',
      'ENSG00000000666': '66',
    }

    self.assertEqual(read_id_map(options), expected)

  def test_symbols(self):
    options = self.get_test_options('symbol')
    expected = {
      'AAA': '11',
      'BBB': '22',
      'CCC': '33',
      'EEE': '55',
      'FFF': '66',
      'aaa': '11',
      'a1': '11',
      'bb': '22',
      'aa': '33',
      'ee': '55',
    }

    self.assertEqual(read_id_map(options), expected)

class MapIdentifiers(unittest.TestCase):
  def test_id_mapping_only(self):
    ids = ['AAA', 'AAB', 'bbb']
    id_map = {
      'AAA': '11',
      'BBB': '22',
      'CCC': '33',
      'EEE': '55',
      'FFF': '66',
    }
    interactions = {}

    expected_mapped_ids = {
      '11': 'AAA',
      '22': 'bbb',
    }
    expected_mapped_interactions = {}

    actual_mapped_ids, actiual_mapped_interactions = map_identifiers(ids, interactions, id_map)
    self.assertEqual(actual_mapped_ids, expected_mapped_ids)
    self.assertEqual(actiual_mapped_interactions, expected_mapped_interactions)

  def test_id_and_interaction_mapping(self):
    ids = ['AAA', 'AAB', 'bbb']
    id_map = {
      'AAA': '11',
      'BBB': '22',
      'CCC': '33',
      'EEE': '55',
      'FFF': '66',
    }
    interactions = {
      'AAA': ['BBB', 'CCC', 'DDD'],
      'AAB': ['BBB', 'CCC', 'EEE'],
      'bbb': ['AAA', 'EEE'],
    }

    expected_mapped_ids = {
      '11': 'AAA',
      '22': 'bbb',
    }
    expected_mapped_interactions = {
      ('11', 'AAA'): {
        '22': { 'symbol': 'BBB' },
        '33': { 'symbol': 'CCC' },
      },
      ('22', 'bbb'): {
        '11': { 'symbol': 'AAA' },
        '55': { 'symbol': 'EEE' },
      }
    }

    actual_mapped_ids, actiual_mapped_interactions = map_identifiers(ids, interactions, id_map)
    self.assertEqual(actual_mapped_ids, expected_mapped_ids)
    self.assertEqual(actiual_mapped_interactions, expected_mapped_interactions)

class CreatePostData(unittest.TestCase):
  def test(self):
    ids = ['111', '333']

    class Options:
      access_key = 'abc123'
      evidence_list = 'IP|APMS'
      include_evidence = True
      include_primary_interactions = True
      include_secondary_interactions = False
      inter_species_excluded = True
      max = 500
      self_interactions_excluded = True
      tax_id = '9606|10090'
      throughput_tag = 'high'
    options = Options()

    expected = {
      'accessKey': 'abc123',
      'additionalIdentifierTypes': 'ENTREZ_GENE',
      'evidenceList': 'IP|APMS',
      'format': 'json',
      'geneList': '111|333',
      'includeEvidence': True,
      'includeInteractors': True,
      'includeInteractorInteractions': False,
      'interSpeciesExcluded': True,
      'max': 500,
      'selfInteractionsExcluded': True,
      'taxId': '9606|10090',
      'throughputTag': 'high',
    }

    self.assertEqual(create_post_data(options, ids), expected)

class WriteBiogridData(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()

  def test(self):
    interactions = [
      { 'BIOGRID_INTERACTION_ID': 1, 'ENTREZ_GENE_A': '111', 'ENTREZ_GENE_B': '222' },
      { 'BIOGRID_INTERACTION_ID': 2, 'ENTREZ_GENE_A': '333', 'ENTREZ_GENE_B': '444' },
      { 'BIOGRID_INTERACTION_ID': 3, 'ENTREZ_GENE_A': '555', 'ENTREZ_GENE_B': '666' },
    ]

    expected = (
      'BIOGRID_INTERACTION_ID\tENTREZ_GENE_A\tENTREZ_GENE_B\n'
      '1\t111\t222\n'
      '2\t333\t444\n'
      '3\t555\t666\n'
    )

    write_biogrid_data(interactions)
    with open('./biogrid-interactions.txt', 'r') as f:
      actual = f.read()
    self.assertEqual(actual, expected)

class ExtractInterationPair(unittest.TestCase):
  def test(self):
    ids = ['111', '222']
    biogrid_data = [
      { 'ENTREZ_GENE_A': '111', 'ENTREZ_GENE_B': '222', 'OFFICIAL_SYMBOL_A': 'AAA', 'OFFICIAL_SYMBOL_B': 'BBB' },
      { 'ENTREZ_GENE_A': '111', 'ENTREZ_GENE_B': '333', 'OFFICIAL_SYMBOL_A': 'AAA', 'OFFICIAL_SYMBOL_B': 'CCC' },
      { 'ENTREZ_GENE_A': '555', 'ENTREZ_GENE_B': '111', 'OFFICIAL_SYMBOL_A': 'EEE', 'OFFICIAL_SYMBOL_B': 'AAA' },
      { 'ENTREZ_GENE_A': '222', 'ENTREZ_GENE_B': '111', 'OFFICIAL_SYMBOL_A': 'BBB', 'OFFICIAL_SYMBOL_B': 'AAA' },
      { 'ENTREZ_GENE_A': '666', 'ENTREZ_GENE_B': '111', 'OFFICIAL_SYMBOL_A': 'FFF', 'OFFICIAL_SYMBOL_B': 'AAA' },
      { 'ENTREZ_GENE_A': '222', 'ENTREZ_GENE_B': '555', 'OFFICIAL_SYMBOL_A': 'BBB', 'OFFICIAL_SYMBOL_B': 'EEE' },
      { 'ENTREZ_GENE_A': '555', 'ENTREZ_GENE_B': '666', 'OFFICIAL_SYMBOL_A': 'EEE', 'OFFICIAL_SYMBOL_B': 'FFF' },
    ]

    expected = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB' },
        '333': { 'symbol': 'CCC' },
        '555': { 'symbol': 'EEE' },
        '666': { 'symbol': 'FFF' },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'AAA' },
        '555': { 'symbol': 'EEE' },
      },
      ('555', 'EEE'): {
        '666': { 'symbol': 'FFF' },
      },
    }

    self.assertEqual(extract_interaction_pairs(biogrid_data, ids), expected)

class MergeInputInteractions(unittest.TestCase):
  def test_no_input_interactions(self):
    input_interactions = {}
    interactions = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB' },
        '333': { 'symbol': 'CCC' },
        '555': { 'symbol': 'EEE' },
        '666': { 'symbol': 'FFF' },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'AAA' },
        '555': { 'symbol': 'EEE' },
      },
      ('555', 'EEE'): {
        '666': { 'symbol': 'FFF' },
      },
    }

    expected = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB', 'known': True, 'isprey': False },
        '333': { 'symbol': 'CCC', 'known': True, 'isprey': False },
        '555': { 'symbol': 'EEE', 'known': True, 'isprey': False },
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'AAA', 'known': True, 'isprey': False },
        '555': { 'symbol': 'EEE', 'known': True, 'isprey': False },
      },
      ('555', 'EEE'): {
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
    }
    self.assertEqual(merge_input_interactions(interactions, input_interactions), expected)

  def test_with_input_interactions(self):
    input_interactions = {
      ('111', 'AAA'): {
        '555': { 'symbol': 'ee' },
        '666': { 'symbol': 'FFF' },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'aa' },
        '666': { 'symbol': 'FFF' },
      }
    }
    interactions = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB' },
        '333': { 'symbol': 'CCC' },
        '555': { 'symbol': 'EEE' },
        '666': { 'symbol': 'FFF' },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'AAA' },
        '555': { 'symbol': 'EEE' },
      },
      ('555', 'EEE'): {
        '666': { 'symbol': 'FFF' },
      },
    }

    expected = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB', 'known': True, 'isprey': False },
        '333': { 'symbol': 'CCC', 'known': True, 'isprey': False },
        '555': { 'symbol': 'ee', 'known': True, 'isprey': True, },
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': True, },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'aa', 'known': True, 'isprey': True, },
        '555': { 'symbol': 'EEE', 'known': True, 'isprey': False },
        '666': { 'symbol': 'FFF', 'known': False, 'isprey': True, },
      },
      ('555', 'EEE'): {
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
    }
    self.assertEqual(merge_input_interactions(interactions, input_interactions), expected)

class ConsolidateTargetSymbols(unittest.TestCase):
  def test_no_input_interactions(self):
    id_type = 'symbol'
    input_ids = {
      '111': 'AAA',
      '222': 'BBB',
    }
    input_interactions = {}
    interactions = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB', 'known': True, 'isprey': False },
        '333': { 'symbol': 'CCC', 'known': True, 'isprey': False },
        '555': { 'symbol': 'ee', 'known': True, 'isprey': False },
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'aa', 'known': True, 'isprey': False },
        '555': { 'symbol': 'EEE', 'known': True, 'isprey': False },
      },
      ('555', 'EEE'): {
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
    }

    expected = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB', 'known': True, 'isprey': False },
        '333': { 'symbol': 'CCC', 'known': True, 'isprey': False },
        '555': { 'symbol': 'ee', 'known': True, 'isprey': False },
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'AAA', 'known': True, 'isprey': False },
        '555': { 'symbol': 'EEE', 'known': True, 'isprey': False },
      },
      ('555', 'EEE'): {
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
    }
    self.assertEqual(consolidate_symbols(id_type, interactions, input_ids, input_interactions), expected)

  def test_with_input_interactions(self):
    id_type = 'symbol'
    input_ids = {
      '111': 'AAA',
      '222': 'BBB',
    }
    input_interactions = {
      ('111', 'AAA'): {
        '555': { 'symbol': 'ee' },
        '666': { 'symbol': 'FFF' },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'aa' },
        '666': { 'symbol': 'FFF' },
      }
    }
    interactions = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB', 'known': True, 'isprey': False },
        '333': { 'symbol': 'CCC', 'known': True, 'isprey': False },
        '555': { 'symbol': 'ee', 'known': True, 'isprey': True, },
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': True, },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'aa', 'known': True, 'isprey': True, },
        '555': { 'symbol': 'EEE', 'known': True, 'isprey': False },
        '666': { 'symbol': 'FFF', 'known': False, 'isprey': True, },
      },
      ('555', 'EEE'): {
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
    }

    expected = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB', 'known': True, 'isprey': False },
        '333': { 'symbol': 'CCC', 'known': True, 'isprey': False },
        '555': { 'symbol': 'ee', 'known': True, 'isprey': True, },
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': True, },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'AAA', 'known': True, 'isprey': True, },
        '555': { 'symbol': 'ee', 'known': True, 'isprey': False },
        '666': { 'symbol': 'FFF', 'known': False, 'isprey': True, },
      },
      ('555', 'ee'): {
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
    }
    self.assertEqual(consolidate_symbols(id_type, interactions, input_ids, input_interactions), expected)

  def test_id_type_not_symbol(self):
    id_type = 'refseqp'
    input_ids = {
      '111': 'NP_11111',
      '222': 'NP_22222',
    }
    input_interactions = {}
    interactions = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB', 'known': True, 'isprey': False },
        '333': { 'symbol': 'CCC', 'known': True, 'isprey': False },
        '555': { 'symbol': 'ee', 'known': True, 'isprey': False },
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'aa', 'known': True, 'isprey': False },
        '555': { 'symbol': 'EEE', 'known': True, 'isprey': False },
      },
      ('555', 'EEE'): {
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
    }

    expected = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB', 'known': True, 'isprey': False },
        '333': { 'symbol': 'CCC', 'known': True, 'isprey': False },
        '555': { 'symbol': 'ee', 'known': True, 'isprey': False },
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'aa', 'known': True, 'isprey': False },
        '555': { 'symbol': 'EEE', 'known': True, 'isprey': False },
      },
      ('555', 'EEE'): {
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
    }
    self.assertEqual(consolidate_symbols(id_type, interactions, input_ids, input_interactions), expected)

class WriteInteractions(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()

  def test_input_saint(self):
    input_ids = {
      '111': 'AAA',
      '222': 'BBB',
    }
    interactions = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB', 'known': True, 'isprey': False },
        '333': { 'symbol': 'CCC', 'known': True, 'isprey': False },
        '555': { 'symbol': 'ee', 'known': True, 'isprey': True, },
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': True, },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'AAA', 'known': True, 'isprey': True, },
        '555': { 'symbol': 'ee', 'known': True, 'isprey': False },
        '666': { 'symbol': 'FFF', 'known': False, 'isprey': True, },
      },
      ('555', 'ee'): {
        '666': { 'symbol': 'FFF', 'known': True, 'isprey': False },
      },
    }
    class Options:
      id_type = 'symbol'
      include_saint_interactions = True
      is_saint = True
    options = Options()

    expected = (
      'source\tsource Entrez\ttarget\ttarget Entrez\tis target a prey\tis target known\n'
      'AAA\t111\tBBB\t222\tFalse\tTrue\n'
      'AAA\t111\tCCC\t333\tFalse\tTrue\n'
      'AAA\t111\tee\t555\tTrue\tTrue\n'
      'AAA\t111\tFFF\t666\tTrue\tTrue\n'
      'BBB\t222\tAAA\t111\tTrue\tTrue\n'
      'BBB\t222\tee\t555\tFalse\tTrue\n'
      'BBB\t222\tFFF\t666\tTrue\tFalse\n'
      'ee\t555\tFFF\t666\tFalse\tTrue\n'
    )

    write_interactions(interactions, input_ids, options)
    with open('./cytoscape.txt', 'r') as f:
      actual = f.read()
    self.assertEqual(actual, expected)

  def test_input_symbol(self):
    input_ids = {
      '111': 'AAA',
      '222': 'BBB',
    }
    interactions = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB' },
        '333': { 'symbol': 'CCC' },
        '555': { 'symbol': 'ee' },
        '666': { 'symbol': 'FFF' },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'AAA' },
        '555': { 'symbol': 'ee' },
        '666': { 'symbol': 'FFF' },
      },
      ('555', 'ee'): {
        '666': { 'symbol': 'FFF' },
      },
    }
    class Options:
      id_type = 'symbol'
      include_saint_interactions = False
      is_saint = False
    options = Options()

    expected = (
      'source\tsource Entrez\ttarget\ttarget Entrez\n'
      'AAA\t111\tBBB\t222\n'
      'AAA\t111\tCCC\t333\n'
      'AAA\t111\tee\t555\n'
      'AAA\t111\tFFF\t666\n'
      'BBB\t222\tAAA\t111\n'
      'BBB\t222\tee\t555\n'
      'BBB\t222\tFFF\t666\n'
      'ee\t555\tFFF\t666\n'
    )

    write_interactions(interactions, input_ids, options)
    with open('./cytoscape.txt', 'r') as f:
      actual = f.read()
    self.assertEqual(actual, expected)

  def test_input_non_symbol(self):
    input_ids = {
      '111': 'NP_11111',
      '222': 'NP_22222',
    }
    interactions = {
      ('111', 'AAA'): {
        '222': { 'symbol': 'BBB' },
        '333': { 'symbol': 'CCC' },
        '555': { 'symbol': 'ee' },
        '666': { 'symbol': 'FFF' },
      },
      ('222', 'BBB'): {
        '111': { 'symbol': 'AAA' },
        '555': { 'symbol': 'ee' },
        '666': { 'symbol': 'FFF' },
      },
      ('555', 'ee'): {
        '666': { 'symbol': 'FFF' },
      },
    }
    class Options:
      id_type = 'refseqp'
      include_saint_interactions = False
      is_saint = False
    options = Options()

    expected = (
      'source ID\tsource symbol\tsource Entrez\ttarget symbol\ttarget Entrez\n'
      'NP_11111\tAAA\t111\tBBB\t222\n'
      'NP_11111\tAAA\t111\tCCC\t333\n'
      'NP_11111\tAAA\t111\tee\t555\n'
      'NP_11111\tAAA\t111\tFFF\t666\n'
      'NP_22222\tBBB\t222\tAAA\t111\n'
      'NP_22222\tBBB\t222\tee\t555\n'
      'NP_22222\tBBB\t222\tFFF\t666\n'
      '\tee\t555\tFFF\t666\n'
    )

    write_interactions(interactions, input_ids, options)
    with open('./cytoscape.txt', 'r') as f:
      actual = f.read()
    self.assertEqual(actual, expected)
 