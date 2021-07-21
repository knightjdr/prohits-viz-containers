import pandas as pd
import pandas.testing as pd_testing
import pyfakefs.fake_filesystem_unittest
import unittest

from .enrich import (
  bh_correction,
  calculate_enrichment,
  count_domains_by_bait,
  filter_saint,
  fishers_test,
  get_background,
  map_file_ids,
  parse_domains,
  read_domains,
  read_gene_map,
  read_saint,
)

class ReadDomains(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()

  def test(self):
    file_contents = (
      '{\n'
      '\t"1": [\n'
      '\t\t{"name": "dA", "start": 10, "end": 21},\n'
      '\t\t{"name": "dB", "start": 30, "end": 45}\n'
      '\t],\n'
      '\t"2": [\n'
      '\t\t{"name": "dC", "start": 15, "end": 35}\n'
      '\t],\n'
      '\t"3": [\n'
      '\t\t{"name": "dA", "start": 20, "end": 31},\n'
      '\t\t{"name": "dD", "start": 35, "end": 65}\n'
      '\t],\n'
      '\t"5": [\n'
      '\t\t{"name": "dD", "start": 10, "end": 40},\n'
      '\t\t{"name": "dE", "start": 50, "end": 55}\n'
      '\t],\n'
      '\t"6": [\n'
      '\t\t{"name": "dF", "start": 30, "end": 45},\n'
      '\t\t{"name": "dA", "start": 60, "end": 71}\n'
      '\t]\n'
      '}\n'
    )
    filepath = '/test/domains.json'
    self.fs.create_file(filepath, contents=file_contents)
  
    expected = {
      '1': [
        { 'name': 'dA', 'start': 10, 'end': 21 },
        { 'name': 'dB', 'start': 30, 'end': 45 },
      ],
      '2': [
        { 'name': 'dC', 'start': 15, 'end': 35 },
      ],
      '3': [
        { 'name': 'dA', 'start': 20, 'end': 31 },
        { 'name': 'dD', 'start': 35, 'end': 65 },
      ],
      '5': [
        { 'name': 'dD', 'start': 10, 'end': 40 },
        { 'name': 'dE', 'start': 50, 'end': 55 },
      ],
      '6': [
        { 'name': 'dF', 'start': 30, 'end': 45 },
        { 'name': 'dA', 'start': 60, 'end': 71 },
      ],
    }

    self.assertEqual(read_domains(filepath), expected)

class ReadGeneMap(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()

  def get_test_options(self, arg_idtype):
    file_contents = (
      '{\n'
      '"1": {\n'
      '"entrez": "11",\n'
      '"refseqp": ["NP_11111"]\n'
      '},\n'
      '"2": {\n'
      '"entrez": "22",\n'
      '"refseqp": ["NP_22222", "NP_02222"]\n'
      '},\n'
      '"3": {\n'
      '"entrez": "33",\n'
      '"refseqp": ["NP_33333"]\n'
      '},\n'
      '"5": {\n'
      '"entrez": "55",\n'
      '"refseqp": ["NP_55555"]\n'
      '},\n'
      '"6": {\n'
      '"entrez": "66",\n'
      '"refseqp": ["NP_66666"]\n'
      '}\n'
      '}\n'
    )
    filepath = '/test/genemap.json'
    self.fs.create_file(filepath, contents=file_contents)

    class Options:
      genemap = filepath
      idtype = arg_idtype

    return Options()

  def test_list_ids(self):
    options = self.get_test_options('refseqp')
    expected = {
      'NP_11111': '1',
      'NP_22222': '2',
      'NP_02222': '2',
      'NP_33333': '3',
      'NP_55555': '5',
      'NP_66666': '6',
    }

    self.assertEqual(read_gene_map(options), expected)

  def test_string_id(self):
    options = self.get_test_options('entrez')
    expected = {
      '11': '1',
      '22': '2',
      '33': '3',
      '55': '5',
      '66': '6',
    }

    self.assertEqual(read_gene_map(options), expected)

class ReadSaint(pyfakefs.fake_filesystem_unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e

  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
    self.setUpPyfakefs()

  def test(self):
    file_contents = (
      'Bait\tPrey\tPreyGene\tSpec\tAvgSpec\tBFDR\n'
      'AAA\tNP_11111\tprey1\t\t10\t0.01\n'
      'AAA\tNP_22222\tprey2\t\t20\t0\n'
      'AAA\tNP_33333\tprey3\t\t30\t0.02\n'
      'AAA\tNP_44444\tprey4\t\t15\t0.01\n'
      'AAA\tNP_55555\tprey5\t\t25\t0.01\n'
      'AAA\tNP_66666\tprey6\t\t40\t0.01\n'
      'BBB\tNP_11111\tprey1\t\t10\t0.05\n'
      'BBB\tNP_22222\tprey2\t\t20\t0.01\n'
      'BBB\tNP_77777\tprey7\t\t30\t0.01\n'
    )
    filepath = '/test/saint.txt'
    self.fs.create_file(filepath, contents=file_contents)
  
    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'NP_11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'NP_33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'BFDR': 0.02 },
      { 'Bait': 'AAA', 'Prey': 'NP_44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.05 },
      { 'Bait': 'BBB', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])

    self.assertEqual(read_saint(filepath), expected)

class MapFileIds(pyfakefs.fake_filesystem_unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e

  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

  def test(self):
    df = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'NP_11111.1', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'NP_33333.3', 'PreyGene': 'prey3', 'AvgSpec': 30, 'BFDR': 0.02 },
      { 'Bait': 'AAA', 'Prey': 'NP_44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.05 },
      { 'Bait': 'BBB', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])
    genemap = {
      'NP_11111': '1',
      'NP_22222': '2',
      'NP_02222': '2',
      'NP_33333': '3',
      'NP_55555': '5',
      'NP_66666': '6',
    }

    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': '1', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': '2', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': '3', 'PreyGene': 'prey3', 'AvgSpec': 30, 'BFDR': 0.02 },
      { 'Bait': 'AAA', 'Prey': 'NP_44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': '5', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': '6', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': '1', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.05 },
      { 'Bait': 'BBB', 'Prey': '2', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])

    self.assertEqual(map_file_ids(df, genemap), expected)

class FilterSaint(pyfakefs.fake_filesystem_unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e

  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

  def get_test_options(self, arg_top_preys):
    df = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'NP_11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'NP_33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'BFDR': 0.02 },
      { 'Bait': 'AAA', 'Prey': 'NP_44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.05 },
      { 'Bait': 'BBB', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])

    class Options:
      fdr = 0.01
      top_preys = arg_top_preys

    return Options(), df

  def test(self):
    options, saint = self.get_test_options(0)
    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'NP_11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'NP_44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])

    self.assertEqual(filter_saint(options, saint), expected)

  def test_top_preys(self):
    options, saint = self.get_test_options(4)
    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'NP_66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'NP_44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
    ])

    self.assertEqual(filter_saint(options, saint), expected)

class GetBackground(pyfakefs.fake_filesystem_unittest.TestCase):
  def get_test_options(self, arg_background):
    domains = {
      '1': [
        { 'name': 'dA', 'start': 10, 'end': 21 },
        { 'name': 'dB', 'start': 30, 'end': 45 },
      ],
      '2': [
        { 'name': 'dC', 'start': 15, 'end': 35 },
      ],
      '3': [
        { 'name': 'dA', 'start': 20, 'end': 31 },
        { 'name': 'dD', 'start': 35, 'end': 65 },
      ],
      '5': [
        { 'name': 'dD', 'start': 10, 'end': 40 },
        { 'name': 'dE', 'start': 50, 'end': 55 },
      ],
      '6': [
        { 'name': 'dF', 'start': 30, 'end': 45 },
        { 'name': 'dA', 'start': 60, 'end': 71 },
      ],
    }
    saint = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'NP_11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'NP_33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'BFDR': 0.02 },
      { 'Bait': 'AAA', 'Prey': 'NP_44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'NP_66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.05 },
      { 'Bait': 'BBB', 'Prey': 'NP_22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])

    class Options:
      background = arg_background

    return Options(), saint, domains
  
  def test_background_from_file(self):
    options, saint, domains = self.get_test_options('file')

    expected = ['NP_11111', 'NP_22222', 'NP_33333', 'NP_44444', 'NP_55555', 'NP_66666', 'NP_77777']
    self.assertListEqual(get_background(options, saint, domains), expected) 

  def test_background_from_all_annotated_genes(self):
    options, saint, domains = self.get_test_options('all')

    expected = ['1', '2', '3', '5', '6']
    self.assertListEqual(get_background(options, saint, domains), expected) 

class ReadSequenceElements(pyfakefs.fake_filesystem_unittest.TestCase):
  def test(self):
    domains = {
      '1': [
        { 'name': 'dA', 'start': 10, 'end': 21 },
        { 'name': 'dB', 'start': 30, 'end': 45 },
        { 'name': 'dA', 'start': 50, 'end': 61 },
      ],
      '2': [
        { 'name': 'dC', 'start': 15, 'end': 35 },
      ],
      '3': [
        { 'name': 'dA', 'start': 20, 'end': 31 },
        { 'name': 'dD', 'start': 35, 'end': 65 },
      ],
      '5': [
        { 'name': 'dD', 'start': 10, 'end': 40 },
        { 'name': 'dE', 'start': 50, 'end': 55 },
      ],
      '6': [
        { 'name': 'dF', 'start': 30, 'end': 45 },
        { 'name': 'dA', 'start': 60, 'end': 71 },
      ],
    }

    background = ['1', '2', '3', '5']

    expected_ids_by_domain = {
      'dA': ['1', '3'],
      'dB': ['1'],
      'dC': ['2'],
      'dD': ['3', '5'],
      'dE': ['5'],
    }
    expected_domains_by_id = {
      '1': {
        'dA': { 'count': 2, 'length': 24 },
        'dB': { 'count': 1, 'length': 16 },
      },
      '2': {
        'dC': { 'count': 1, 'length': 21 },
      },
      '3': {
        'dA': { 'count': 1, 'length': 12 },
        'dD': { 'count': 1, 'length': 31 },
      },
      '5': {
        'dD': { 'count': 1, 'length': 31 },
        'dE': { 'count': 1, 'length': 6 },
      },
    }

    acutal_domains_by_id, acutal_ids_by_domain = parse_domains(domains, background)

    self.assertEqual(acutal_ids_by_domain, expected_ids_by_domain)
    self.assertEqual(acutal_domains_by_id, expected_domains_by_id)

class CountDomainsByBait(unittest.TestCase):
  def test(self):
    saint = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': '1', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': '2', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': '3', 'PreyGene': 'prey3', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': '5', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': '6', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': '2', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'NP_77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])
    domains_by_id = {
      '1': {
        'dA': { 'count': 2, 'length': 24 },
        'dB': { 'count': 1, 'length': 16 },
      },
      '2': {
        'dC': { 'count': 1, 'length': 21 },
      },
      '3': {
        'dA': { 'count': 1, 'length': 12 },
        'dD': { 'count': 1, 'length': 31 },
      },
      '5': {
        'dD': { 'count': 1, 'length': 31 },
        'dE': { 'count': 1, 'length': 6 },
      },
    }

    expected = {
      'AAA': {
        'preysInDatabase': 4,
        'domain': {
          'dA': {
            'countByAccession': [2, 1],
            'lengthByAccession': [24, 12],
            'preys': ['prey1', 'prey3'],
          },
          'dB': {
            'countByAccession': [1],
            'lengthByAccession': [16],
            'preys': ['prey1'],
          },
          'dC': {
            'countByAccession': [1],
            'lengthByAccession': [21],
            'preys': ['prey2'],
          },
          'dD': {
            'countByAccession': [1, 1],
            'lengthByAccession': [31, 31],
            'preys': ['prey3', 'prey5'],
          },
          'dE': {
            'countByAccession': [1],
            'lengthByAccession': [6],
            'preys': ['prey5'],
          },
        },
      },
      'BBB': {
        'preysInDatabase': 1,
        'domain': {
          'dC': {
            'countByAccession': [1],
            'lengthByAccession': [21],
            'preys': ['prey2'],
          },
        }
      },
    }

    self.assertEqual(count_domains_by_bait(saint, domains_by_id), expected)

class FishersTest(unittest.TestCase):
  def test(self):
    n11 = 10
    n1p = 15
    np1 = 15
    npp = 100

    expected = 4e-7

    self.assertAlmostEqual(fishers_test(n11, n1p, np1, npp), expected, places=5)

class BhCorrection(unittest.TestCase):
  def test(self):
    fdr = 0.01
    pvalues = {
      'd1': 0.01,
      'd2': 0.00001,
      'd3': 0.001,
      'd4': 0.0001,
      'd5': 0.001,
    }

    expected_adj_pvalues = {
      'd1': 0.0125000,
      'd5': 0.0016666666666666668,
      'd3': 0.0016666666666666668,
      'd4': 0.0002500,
      'd2': 0.0000500,
    }
    expected_corrected_fdr = {
      'd1': 0.008,
      'd5': 0.006,
      'd3': 0.006,
      'd4': 0.004,
      'd2': 0.002,
    }

    actual_adj_pvalues, actual_corrected_fdr = bh_correction(pvalues, fdr)
    self.assertEqual(actual_adj_pvalues, expected_adj_pvalues)
    self.assertEqual(actual_corrected_fdr, expected_corrected_fdr)

class CalculateEnrichment(unittest.TestCase):
  def test(self):
    background_size = 1000
    fdr = 0.01
    domains_by_bait = {
      'AAA': {
        'preysInDatabase': 3,
        'domain': {
          'dA': {
            'countByAccession': [1, 2],
            'lengthByAccession': [12, 31],
            'preys': ['prey1', 'prey3'],
          },
          'dB': {
            'countByAccession': [1],
            'lengthByAccession': [11],
            'preys': ['prey3'],
          },
        },
      },
      'BBB': {
        'preysInDatabase': 2,
        'domain': {
          'dA': {
            'countByAccession': [2],
            'lengthByAccession': [31],
            'preys': ['prey3'],
          },
          'dB': {
            'countByAccession': [1],
            'lengthByAccession': [11],
            'preys': ['prey3'],
          },
        },
      },
    }
    ids_by_domain = {
      'dA': ['1', '3', '4', '5'],
      'dB': ['3'],
    }
  
    expected = {
      'AAA': [
        {
          'domain': 'dA',
          'no_genes_with_domain': 2,
          'no_genes': 3,
          'background_size_w_domain': 4,
          'background_size': 1000,
          'fold_enrichment': 166.66666666666666,
          'pvalue': 3.5987891699293464e-05,
          'adj_pvalue': 7.197578339858693e-05,
          'bh_fdr': 0.005,
          'genes': ['prey1', 'prey3'],
        },
        {
          'domain': 'dB',
          'no_genes_with_domain': 1,
          'no_genes': 3,
          'background_size_w_domain': 1,
          'background_size': 1000,
          'fold_enrichment': 333.3333333333333,
          'pvalue': 0.0029999999999970202,
          'adj_pvalue': 0.0029999999999970202,
          'bh_fdr': 0.01,
          'genes': ['prey3'],
        },
      ],
      'BBB': [
        {
          'domain': 'dB',
          'no_genes_with_domain': 1,
          'no_genes': 2,
          'background_size_w_domain': 1,
          'background_size': 1000,
          'fold_enrichment': 500.0,
          'pvalue': 0.0019999999999980147,
          'adj_pvalue': 0.003999999999996029,
          'bh_fdr': 0.005,
          'genes': ['prey3'],
        },
        {
          'domain': 'dA',
          'no_genes_with_domain': 1,
          'no_genes': 2,
          'background_size_w_domain': 4,
          'background_size': 1000,
          'fold_enrichment': 125.0,
          'pvalue': 0.007987987987994435,
          'adj_pvalue': 0.007987987987994435,
          'bh_fdr': 0.01,
          'genes': ['prey3'],
        },
      ],
    }
    
    actual = calculate_enrichment(domains_by_bait, ids_by_domain, background_size, fdr)
    self.assertEqual(actual, expected)
    