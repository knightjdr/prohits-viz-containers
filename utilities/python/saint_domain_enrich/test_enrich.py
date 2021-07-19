import pandas as pd
import pandas.testing as pd_testing
import pyfakefs.fake_filesystem_unittest
import unittest

from .enrich import (
  # bh_correction,
  # calculate_enrichment,
  # count_elements_by_bait,
  filter_saint,
  # fishers_test,
  get_background,
  read_domains,
  # read_sequence_elements,
  read_saint,
)

class ReadDomains(pyfakefs.fake_filesystem_unittest.TestCase):
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
      '{\n'
      '"1": ["domainA", "domainB"],\n'
      '"2": ["domainC"],\n'
      '"3": ["domainA", "domainD"],\n'
      '"5": ["domainD", "domainE"],\n'
      '"6": ["domainA", "domainF"]\n'
      '}\n'
    )
    filepath = '/test/domains.json'
    self.fs.create_file(filepath, contents=file_contents)
  
    expected = {
      '1': ['domainA', 'domainB'],
      '2': ['domainC'],
      '3': ['domainA', 'domainD'],
      '5': ['domainD', 'domainE'],
      '6': ['domainA', 'domainF'],
    }

    self.assertEqual(read_domains(filepath), expected)

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
      'AAA\tP11111\tprey1\t\t10\t0.01\n'
      'AAA\tP22222\tprey2\t\t20\t0\n'
      'AAA\tP33333\tprey3\t\t30\t0.02\n'
      'AAA\tP44444\tprey4\t\t15\t0.01\n'
      'AAA\tP55555\tprey5\t\t25\t0.01\n'
      'AAA\tP66666\tprey6\t\t40\t0.01\n'
      'BBB\tP11111\tprey1\t\t10\t0.05\n'
      'BBB\tP22222\tprey2\t\t20\t0.01\n'
      'BBB\tP77777\tprey7\t\t30\t0.01\n'
    )
    filepath = '/test/saint.txt'
    self.fs.create_file(filepath, contents=file_contents)
  
    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'BFDR': 0.02 },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.05 },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])

    self.assertEqual(read_saint(filepath), expected)

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
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'BFDR': 0.02 },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.05 },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])

    class Options:
      fdr = 0.01
      top_preys = arg_top_preys

    return Options(), df

  def test(self):
    options, saint = self.get_test_options(0)
    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])

    self.assertEqual(filter_saint(options, saint), expected)

  def test_top_preys(self):
    options, saint = self.get_test_options(4)
    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
    ])

    self.assertEqual(filter_saint(options, saint), expected)

class GetBackground(pyfakefs.fake_filesystem_unittest.TestCase):
  def get_test_options(self, arg_background):
    domains = {
      '1': ['domainA', 'domainB'],
      '2': ['domainC'],
      '3': ['domainA', 'domainD'],
      '5': ['domainD', 'domainE'],
      '6': ['domainA', 'domainF'],
    }
    saint = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'BFDR': 0.02 },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.05 },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P77777', 'PreyGene': 'prey7', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])

    class Options:
      background = arg_background

    return Options(), saint, domains
  
  def test_background_from_file(self):
    options, saint, domains = self.get_test_options('file')

    expected = ['P11111', 'P22222', 'P33333', 'P44444', 'P55555', 'P66666', 'P77777']
    self.assertListEqual(get_background(options, saint, domains), expected) 

  def test_background_from_all_annotated_genes(self):
    options, saint, domains = self.get_test_options('all')

    expected = ['1', '2', '3', '5', '6']
    self.assertListEqual(get_background(options, saint, domains), expected) 

""" class BhCorrection(unittest.TestCase):
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
    self.assertEqual(actual_corrected_fdr, expected_corrected_fdr) """

""" class CalculateEnrichment(unittest.TestCase):
  def test(self):
    background_size = 1000
    fdr = 0.01
    elements_by_bait = {
      'AAA': {
        'preysInDatabase': 3,
        'domain': {
          'd-a': {
            'countByAccession': [1, 2],
            'lengthByAccession': [12, 31],
            'preys': ['prey1', 'prey3'],
          },
          'd-b': {
            'countByAccession': [1],
            'lengthByAccession': [11],
            'preys': ['prey3'],
          },
        },
        'region': {},
      },
      'BBB': {
        'preysInDatabase': 2,
        'domain': {
          'd-a': {
            'countByAccession': [2],
            'lengthByAccession': [31],
            'preys': ['prey3'],
          },
          'd-b': {
            'countByAccession': [1],
            'lengthByAccession': [11],
            'preys': ['prey3'],
          },
        },
        'region': {},
      },
    }
    sequence_elements = {
      'd-a': { 'P11111': 'prey1', 'P33333': 'prey3', 'P44444': 'prey4', 'P55555': 'prey5' },
      'd-b': { 'P33333': 'prey3' },
    }
  
    expected = {
      'AAA': {
        'domain': [
          {
            'element': 'd-a',
            'no_genes_with_element': 2,
            'no_genes': 3,
            'background_size_w_element': 4,
            'background_size': 1000,
            'fold_enrichment': 166.66666666666666,
            'pvalue': 3.5987891699293464e-05,
            'adj_pvalue': 7.197578339858693e-05,
            'bh_fdr': 0.005,
            'genes': ['prey1', 'prey3'],
          },
          {
            'element': 'd-b',
            'no_genes_with_element': 1,
            'no_genes': 3,
            'background_size_w_element': 1,
            'background_size': 1000,
            'fold_enrichment': 333.3333333333333,
            'pvalue': 0.0029999999999970202,
            'adj_pvalue': 0.0029999999999970202,
            'bh_fdr': 0.01,
            'genes': ['prey3'],
          },
        ],
        'region': [],
      },
      'BBB': {
        'domain': [
          {
            'element': 'd-b',
            'no_genes_with_element': 1,
            'no_genes': 2,
            'background_size_w_element': 1,
            'background_size': 1000,
            'fold_enrichment': 500.0,
            'pvalue': 0.0019999999999980147,
            'adj_pvalue': 0.003999999999996029,
            'bh_fdr': 0.005,
            'genes': ['prey3'],
          },
          {
            'element': 'd-a',
            'no_genes_with_element': 1,
            'no_genes': 2,
            'background_size_w_element': 4,
            'background_size': 1000,
            'fold_enrichment': 125.0,
            'pvalue': 0.007987987987994435,
            'adj_pvalue': 0.007987987987994435,
            'bh_fdr': 0.01,
            'genes': ['prey3'],
          },
        ],
        'region': [],
      },
    }
    
    actual = calculate_enrichment(elements_by_bait, sequence_elements, background_size, fdr)
    self.assertEqual(actual, expected) """

""" class CountElementsByBait(unittest.TestCase):
  def test(self):
    saint = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1' },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2' },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3' },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2' },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3' },
      { 'Bait': 'BBB', 'Prey': 'P44444', 'PreyGene': 'prey4' },
    ])
    sequence_elements_by_accession = {
      'P11111': {
        'd-a': { 'count': 1, 'length': 12, 'type': 'domain' },
      },
      'P22222': {
        'r-a': { 'count': 1, 'length': 12, 'type': 'region' }
      },
      'P33333': {
        'd-a': { 'count': 2, 'length': 31, 'type': 'domain' },
        'd-b': { 'count': 1, 'length': 11, 'type': 'domain' },
        'r-a': { 'count': 2, 'length': 20, 'type': 'region' },
        'r-b': { 'count': 1, 'length': 10, 'type': 'region' },
      },
    }

    expected = {
      'AAA': {
        'preysInDatabase': 3,
        'domain': {
          'd-a': {
            'countByAccession': [1, 2],
            'lengthByAccession': [12, 31],
            'preys': ['prey1', 'prey3'],
          },
          'd-b': {
            'countByAccession': [1],
            'lengthByAccession': [11],
            'preys': ['prey3'],
          },
        },
        'region': {
          'r-a': {
            'countByAccession': [1, 2],
            'lengthByAccession': [12, 20],
            'preys': ['prey2', 'prey3'],
          },
          'r-b': {
            'countByAccession': [1],
            'lengthByAccession': [10],
            'preys': ['prey3'],
          },
        },
      },
      'BBB': {
        'preysInDatabase': 2,
        'domain': {
          'd-a': {
            'countByAccession': [2],
            'lengthByAccession': [31],
            'preys': ['prey3'],
          },
          'd-b': {
            'countByAccession': [1],
            'lengthByAccession': [11],
            'preys': ['prey3'],
          },
        },
        'region': {
          'r-a': {
            'countByAccession': [1, 2],
            'lengthByAccession': [12, 20],
            'preys': ['prey2', 'prey3'],
          },
          'r-b': {
            'countByAccession': [1],
            'lengthByAccession': [10],
            'preys': ['prey3'],
          },
        },
      },
    }

    self.assertEqual(count_elements_by_bait(saint, sequence_elements_by_accession), expected) """

""" class FishersTest(unittest.TestCase):
  def test(self):
    n11 = 10
    n1p = 15
    np1 = 15
    npp = 100

    expected = 4e-7

    self.assertAlmostEqual(fishers_test(n11, n1p, np1, npp), expected, places=5) """

""" class ReadSequenceElements(pyfakefs.fake_filesystem_unittest.TestCase):
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
      '[\n'
      '\t{\n'
      '\t\t"gene": "gene1",\n'
      '\t\t"protein-expression": { "cells": { "A-549": { "intensity": 0.01 } } },\n'
      '\t\t"uniprot": ["P11111", "Q11111"],\n'
      '\t\t"domains": [\n'
      '\t\t\t{"name": "d1", "type": "domain", "start": 10, "end": 21 }\n'
      '\t\t]\n'
      '\t},\n'
      '\t{\n'
      '\t\t"gene": "gene2",\n'
      '\t\t"protein-expression": { "cells": { "A-549": { "intensity": 0 } } },\n'
      '\t\t"uniprot": ["P22222"],\n'
      '\t\t"domains": [\n'
      '\t\t\t{"name": "d2", "type": "domain", "start": 5, "end": 10}\n'
      '\t\t]\n'
      '\t},\n'
      '\t{\n'
      '\t\t"gene": "gene3",\n'
      '\t\t"protein-expression": { "cells": { "A-549": { "intensity": 1 } } },\n'
      '\t\t"uniprot": ["P33333"],\n'
      '\t\t"domains": []\n'
      '\t},\n'
      '\t{\n'
      '\t\t"gene": "gene4",\n'
      '\t\t"protein-expression": { "cells": {} },\n'
      '\t\t"uniprot": ["P44444"],\n'
      '\t\t"domains": [\n'
      '\t\t\t{"name": "d4", "type": "domain", "start": 10, "end": 21},\n'
      '\t\t\t{"name": "d4-b", "type": "domain", "start": 10, "end": 20},\n'
      '\t\t\t{"name": "r4", "type": "region", "start": 10, "end": 25},\n'
      '\t\t\t{"name": "d4", "type": "domain", "start": 12, "end": 30},\n'
      '\t\t\t{"name": "r4-b", "type": "region", "start": 12, "end": 21},\n'
      '\t\t\t{"name": "r4", "type": "region", "start": 12, "end": 15}\n'
      '\t\t]\n'
      '\t}\n'
      ']\n'
    )
    filepath = '/test/gix.json'
    self.fs.create_file(filepath, contents=file_contents)

    background = ['Q11111', 'P33333', 'P44444']

    expected_sequence_elements = {
      'd1': { 'P11111': 'gene1' },
      'd4': { 'P44444': 'gene4' },
      'd4-b': { 'P44444': 'gene4' },
      'r4': { 'P44444': 'gene4' },
      'r4-b': { 'P44444': 'gene4' },
    }
    expected_sequence_elements_by_accession = {
      'P11111': {
        'd1': { 'count': 1, 'length': 12, 'type': 'domain' },
      },
      'Q11111': {
        'd1': { 'count': 1, 'length': 12, 'type': 'domain' },
      },
      'P33333': {},
      'P44444': {
        'd4': { 'count': 2, 'length': 31, 'type': 'domain' },
        'd4-b': { 'count': 1, 'length': 11, 'type': 'domain' },
        'r4': { 'count': 2, 'length': 20, 'type': 'region' },
        'r4-b': { 'count': 1, 'length': 10, 'type': 'region' },
      },
    }

    acutal_sequence_elements_by_accesion, acutal_sequence_elements = read_sequence_elements(filepath, background)

    self.assertEqual(acutal_sequence_elements, expected_sequence_elements)
    self.assertEqual(acutal_sequence_elements_by_accesion, expected_sequence_elements_by_accession) """
    