import math
import pandas as pd
import pandas.testing as pd_testing
import pyfakefs.fake_filesystem_unittest
import unittest

from .main import (
  add_specificity_to_saint,
  get_specificty_calculator,
  read_saint,
)

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
      'Bait\tPrey\tPreyGene\tAvgSpec\tSpec\tctrlCounts\n'
      'AAA\tP11111\tprey1\t10\t10|10\t0|0\n'
      'AAA\tP22222\tprey2\t20\t20|20\t5|4\n'
      'AAA\tP33333\tprey3\t30\t30|30\t0|3\n'
      'AAA\tP44444\tprey4\t15\t15|15\t7|8\n'
      'AAA\tP55555\tprey5\t25\t25|25\t0|0\n'
      'AAA\tP66666\tprey6\t40\t40|40\t1|1\n'
      'BBB\tP11111\tprey1\t10\t10|10\t0|0\n'
      'BBB\tP22222\tprey2\t20\t20|20\t5|4\n'
      'BBB\tP33333\tprey3\t30\t30|30\t0|3\n'
    )
    filepath = '/test/saint.txt'
    self.fs.create_file(filepath, contents=file_contents)

    control_subtract = False

    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'Spec': '10|10', 'ctrlCounts': '0|0', 'Abundance': 10, 'Replicates': '10|10' },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'Spec': '20|20', 'ctrlCounts': '5|4', 'Abundance': 20, 'Replicates': '20|20' },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'Spec': '30|30', 'ctrlCounts': '0|3', 'Abundance': 30, 'Replicates': '30|30' },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'Spec': '15|15', 'ctrlCounts': '7|8', 'Abundance': 15, 'Replicates': '15|15' },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'Spec': '25|25', 'ctrlCounts': '0|0', 'Abundance': 25, 'Replicates': '25|25' },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'Spec': '40|40', 'ctrlCounts': '1|1', 'Abundance': 40, 'Replicates': '40|40' },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'Spec': '10|10', 'ctrlCounts': '0|0', 'Abundance': 10, 'Replicates': '10|10' },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'Spec': '20|20', 'ctrlCounts': '5|4', 'Abundance': 20, 'Replicates': '20|20' },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'Spec': '30|30', 'ctrlCounts': '0|3', 'Abundance': 30, 'Replicates': '30|30' },
    ])

    self.assertEqual(read_saint(filepath, control_subtract), expected)

  def test_control_subtract(self):
    file_contents = (
      'Bait\tPrey\tPreyGene\tAvgSpec\tSpec\tctrlCounts\n'
      'AAA\tP11111\tprey1\t10\t10|10\t0|0\n'
      'AAA\tP22222\tprey2\t20\t20|20\t5|4\n'
      'AAA\tP33333\tprey3\t30\t30|30\t0|3\n'
      'AAA\tP44444\tprey4\t15\t15|15\t7|8\n'
      'AAA\tP55555\tprey5\t25\t25|25\t0|0\n'
      'AAA\tP66666\tprey6\t40\t40|40\t1|1\n'
      'BBB\tP11111\tprey1\t10\t10|10\t0|0\n'
      'BBB\tP22222\tprey2\t20\t20|20\t5|4\n'
      'BBB\tP33333\tprey3\t30\t30|30\t0|3\n'
    )
    filepath = '/test/saint.txt'
    self.fs.create_file(filepath, contents=file_contents)

    control_subtract = True

    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'Spec': '10|10', 'ctrlCounts': '0|0', 'Abundance': 10, 'Replicates': '10.0|10.0' },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'Spec': '20|20', 'ctrlCounts': '5|4', 'Abundance': 15.5, 'Replicates': '15.5|15.5' },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'Spec': '30|30', 'ctrlCounts': '0|3', 'Abundance': 28.5, 'Replicates': '28.5|28.5' },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'AvgSpec': 15, 'Spec': '15|15', 'ctrlCounts': '7|8', 'Abundance': 7.5, 'Replicates': '7.5|7.5' },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'AvgSpec': 25, 'Spec': '25|25', 'ctrlCounts': '0|0', 'Abundance': 25, 'Replicates': '25.0|25.0' },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'AvgSpec': 40, 'Spec': '40|40', 'ctrlCounts': '1|1', 'Abundance': 39, 'Replicates': '39.0|39.0' },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'AvgSpec': 10, 'Spec': '10|10', 'ctrlCounts': '0|0', 'Abundance': 10, 'Replicates': '10.0|10.0' },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'AvgSpec': 20, 'Spec': '20|20', 'ctrlCounts': '5|4', 'Abundance': 15.5, 'Replicates': '15.5|15.5' },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3', 'AvgSpec': 30, 'Spec': '30|30', 'ctrlCounts': '0|3', 'Abundance': 28.5, 'Replicates': '28.5|28.5' },
    ])

    self.assertEqual(read_saint(filepath, control_subtract), expected)

class CalculateSpecificity(unittest.TestCase):
  def test_dscore(self):
    metric = 'dscore'
    df = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '7|8' },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'Abundance': 25, 'Replicates': '25|25', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'Abundance': 40, 'Replicates': '40|40', 'ctrlCounts': '1|1' },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'CCC', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|0' },
      { 'Bait': 'CCC', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '5|4' },
      { 'Bait': 'CCC', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|3' },
    ])

    calculate_specificity = get_specificty_calculator(df, metric)

    param_list = [
      ('P11111', 10, '10|10', 3.16),
      ('P33333', 30, '30|30', 5.48),
      ('P44444', 15, '15|15', 11.62),
      ('P22222', 20, '20|20', 4.47),
    ]

    for prey, spec, reps, expected in param_list:
      with self.subTest():
        self.assertEqual(calculate_specificity(prey, spec, reps), expected)

  def test_fc(self):
    metric = 'fe'
    df = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '7|8' },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'Abundance': 25, 'Replicates': '25|25', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'Abundance': 40, 'Replicates': '40|40', 'ctrlCounts': '1|1' },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'CCC', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|0' },
      { 'Bait': 'CCC', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '5|4' },
      { 'Bait': 'CCC', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|3' },
    ])

    calculate_specificity = get_specificty_calculator(df, metric)

    param_list = [
      ('P11111', 10, 0.8),
      ('P33333', 30, 1.33),
      ('P44444', 15, math.inf),
      ('P22222', 20, 1.14),
    ]

    for prey, spec, expected in param_list:
      with self.subTest():
        self.assertEqual(calculate_specificity(prey, spec), expected)

  def test_sscore(self):
    metric = 'sscore'
    df = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '7|8' },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'Abundance': 25, 'Replicates': '25|25', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'Abundance': 40, 'Replicates': '40|40', 'ctrlCounts': '1|1' },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'CCC', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|0' },
      { 'Bait': 'CCC', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '5|4' },
      { 'Bait': 'CCC', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|3' },
    ])

    calculate_specificity = get_specificty_calculator(df, metric)

    param_list = [
      ('P11111', 10, 3.16),
      ('P33333', 30, 5.48),
      ('P44444', 15, 6.71),
      ('P22222', 20, 4.47),
    ]

    for prey, spec, expected in param_list:
      with self.subTest():
        self.assertEqual(calculate_specificity(prey, spec), expected)

  def test_wdscore(self):
    metric = 'wdscore'
    df = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '7|8' },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'Abundance': 25, 'Replicates': '25|25', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'Abundance': 40, 'Replicates': '40|40', 'ctrlCounts': '1|1' },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'CCC', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|0' },
      { 'Bait': 'CCC', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '5|4' },
      { 'Bait': 'CCC', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|3' },
    ])

    calculate_specificity = get_specificty_calculator(df, metric)

    param_list = [
      ('P11111', 10, '10|10', 3.16),
      ('P33333', 30, '30|30', 5.48),
      ('P44444', 15, '15|15', 20.12),
      ('P22222', 20, '20|20', 4.47),
    ]

    for prey, spec, reps, expected in param_list:
      with self.subTest():
        self.assertEqual(calculate_specificity(prey, spec, reps), expected)

  def test_zscore(self):
    metric = 'zscore'
    df = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '7|8' },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'Abundance': 25, 'Replicates': '25|25', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'Abundance': 40, 'Replicates': '40|40', 'ctrlCounts': '1|1' },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'CCC', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|0' },
      { 'Bait': 'CCC', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '5|4' },
      { 'Bait': 'CCC', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|3' },
    ])

    calculate_specificity = get_specificty_calculator(df, metric)

    param_list = [
      ('P11111', 10, -0.58),
      ('P33333', 30, 0.58),
      ('P44444', 15, 1.15),
      ('P22222', 20, 0.58),
    ]

    for prey, spec, expected in param_list:
      with self.subTest():
        self.assertEqual(calculate_specificity(prey, spec), expected)

class AddSpecificityToSaint(unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e

  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

  def test(self):
    metric = 'fc'
    df = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '7|8' },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'Abundance': 25, 'Replicates': '25|25', 'ctrlCounts': '0|0' },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'Abundance': 40, 'Replicates': '40|40', 'ctrlCounts': '1|1' },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0' },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4' },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3' },
      { 'Bait': 'CCC', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|0' },
      { 'Bait': 'CCC', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '5|4' },
      { 'Bait': 'CCC', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|3' },
    ])
    calculate_specificity = get_specificty_calculator(df, metric)

    expected = pd.DataFrame([
      { 'Bait': 'AAA', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0', 'Specificity': 0.8 },
      { 'Bait': 'AAA', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4', 'Specificity': 1.14 },
      { 'Bait': 'AAA', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3', 'Specificity': 1.33 },
      { 'Bait': 'AAA', 'Prey': 'P44444', 'PreyGene': 'prey4', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '7|8', 'Specificity': math.inf },
      { 'Bait': 'AAA', 'Prey': 'P55555', 'PreyGene': 'prey5', 'Abundance': 25, 'Replicates': '25|25', 'ctrlCounts': '0|0', 'Specificity': math.inf },
      { 'Bait': 'AAA', 'Prey': 'P66666', 'PreyGene': 'prey6', 'Abundance': 40, 'Replicates': '40|40', 'ctrlCounts': '1|1', 'Specificity': math.inf },
      { 'Bait': 'BBB', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 10, 'Replicates': '10|10', 'ctrlCounts': '0|0', 'Specificity': 0.8 },
      { 'Bait': 'BBB', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 20, 'Replicates': '20|20', 'ctrlCounts': '5|4', 'Specificity': 1.14 },
      { 'Bait': 'BBB', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 30, 'Replicates': '30|30', 'ctrlCounts': '0|3', 'Specificity': 1.33 },
      { 'Bait': 'CCC', 'Prey': 'P11111', 'PreyGene': 'prey1', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|0', 'Specificity': 1.5 },
      { 'Bait': 'CCC', 'Prey': 'P22222', 'PreyGene': 'prey2', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '5|4', 'Specificity': 0.75 },
      { 'Bait': 'CCC', 'Prey': 'P33333', 'PreyGene': 'prey3', 'Abundance': 15, 'Replicates': '15|15', 'ctrlCounts': '0|3', 'Specificity': 0.5 },
    ])
    self.assertEqual(add_specificity_to_saint(df, calculate_specificity), expected)