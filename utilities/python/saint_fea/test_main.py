import pandas as pd
import pandas.testing as pd_testing
import pyfakefs.fake_filesystem_unittest
import unittest

from .main import (
  create_query_lists,
  read_saint,
)

class CreateQueryLists(unittest.TestCase):
  def test(self):
    df = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1' },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2' },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2' },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3' },
    ])

    expected_query = {
      'AAA': ['P11111', 'P22222'],
      'BBB': ['P22222', 'P33333']
    }
    expected_accessions_to_symbol = {
      'P11111': 'prey1',
      'P22222': 'prey2',
      'P33333': 'prey3',
    }

    actual_query, actual_accessions_to_symbol = create_query_lists(df)
    self.assertEqual(actual_query, expected_query)
    self.assertEqual(actual_accessions_to_symbol, expected_accessions_to_symbol)

class ReadSaint(pyfakefs.fake_filesystem_unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e

  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
    self.setUpPyfakefs()

  def get_test_options(self, arg0):
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
      'BBB\tP33333\tprey3\t\t30\t0.01\n'
    )
    filepath = '/test/saint.txt'
    self.fs.create_file(filepath, contents=file_contents)


    class Options:
      fdr = 0.01
      saint = filepath
      top_preys = arg0

    return Options()

  def test(self):
    options = self.get_test_options(0)
    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'BFDR': 0.01 },
    ])

    self.assertEqual(read_saint(options), expected)

  def test_top_preys(self):
    options = self.get_test_options(4)
    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'BFDR': 0.01 },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0 },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'BFDR': 0.01 },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'BFDR': 0.01 },
    ])

    self.assertEqual(read_saint(options), expected)
  