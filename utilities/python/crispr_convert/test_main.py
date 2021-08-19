import pandas as pd
import pandas.testing as pd_testing
import pyfakefs.fake_filesystem_unittest
import unittest

from .main import (
  convert_mageck,
  convert_ranks,
  extract_condition_from_filename,
  get_column_names,
  get_files,
  merge_and_add_conditions,
  move_condition_column,
)

class ConvertMageck(pyfakefs.fake_filesystem_unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e

  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
    self.setUpPyfakefs()

  def test_mle(self):
    file_contents = (
      'Gene\tsgRNA\tHL60|beta\tHL60|z\tHL60|p-value\tHL60|fdr\tHL60|wald-p-value\tHL60|wald-fdr\tKBM7|beta\tKBM7|z\tKBM7|p-value\tKBM7|fdr\tKBM7|wald-p-value\tKBM7|wald-fdr\n'
      'RNF14\t10\t0.24927\t0.72077\t0.36256\t0.75648\t0.47105\t0.9999\t0.57276\t1.6565\t0.06468\t0.32386\t0.097625\t0.73193\n'
      'RNF10\t10\t0.10159\t0.29373\t0.92087\t0.98235\t0.76896\t0.9999\t0.11341\t0.32794\t0.90145\t0.97365\t0.74296\t0.98421\n'
      'RNF11\t10\t3.6354\t10.513\t0.0002811\t0.021739\t7.5197e-26\t1.3376e-22\t2.5928\t7.4925\t0.0014898\t0.032024\t6.7577e-14\t1.33e-11\n'
    )
    filepath = '/folder/mle.txt'
    self.fs.create_file(filepath, contents=file_contents)

    files = [filepath]

    expected = pd.DataFrame([
      { 'condition': 'HL60', 'Gene': 'RNF14', 'sgRNA': 10, 'beta': 0.24927, 'z': 0.72077, 'p-value': 0.36256, 'fdr': 0.75648, 'wald-p-value': 0.47105, 'wald-fdr': 0.9999 },
      { 'condition': 'HL60', 'Gene': 'RNF10', 'sgRNA': 10, 'beta': 0.10159, 'z': 0.29373, 'p-value': 0.92087, 'fdr': 0.98235, 'wald-p-value': 0.76896, 'wald-fdr': 0.9999 },
      { 'condition': 'HL60', 'Gene': 'RNF11', 'sgRNA': 10, 'beta': 3.6354, 'z': 10.513, 'p-value': 0.0002811, 'fdr': 0.021739, 'wald-p-value': 7.5197e-26, 'wald-fdr': 1.3376e-22 },
      { 'condition': 'KBM7', 'Gene': 'RNF14', 'sgRNA': 10, 'beta': 0.57276, 'z': 1.6565, 'p-value': 0.06468, 'fdr': 0.32386, 'wald-p-value': 0.097625, 'wald-fdr': 0.73193 },
      { 'condition': 'KBM7', 'Gene': 'RNF10', 'sgRNA': 10, 'beta': 0.11341, 'z': 0.32794, 'p-value': 0.90145, 'fdr': 0.97365, 'wald-p-value': 0.74296, 'wald-fdr': 0.98421 },
      { 'condition': 'KBM7', 'Gene': 'RNF11', 'sgRNA': 10, 'beta': 2.5928, 'z': 7.4925, 'p-value': 0.0014898, 'fdr': 0.032024, 'wald-p-value': 6.7577e-14, 'wald-fdr': 1.33e-11 },
    ])
    self.assertEqual(convert_mageck(files), expected)

  def test_rra(self):
    file_content_condition1 = (
      'id	num	neg|score	neg|p-value	neg|fdr	neg|rank	neg|goodsgrna	neg|lfc	pos|score	pos|p-value	pos|fdr	pos|rank	pos|goodsgrna	pos|lfc\n'
      'RPS5	5	5.9353e-10	2.5851e-07	0.000707	2	5	-30.64996	1.0	1.0	1.0	19149	0	0\n'
      'GTF2B	5	2.0462e-10	2.5851e-07	0.000707	1	5	-32.18633	1.0	1.0	1.0	19150	0	0\n'
    )
    filepath_condition1 = '/folder/condition1.txt'
    self.fs.create_file(filepath_condition1, contents=file_content_condition1)
    file_content_condition2 = (
      'id	num	neg|score	neg|p-value	neg|fdr	neg|rank	neg|goodsgrna	neg|lfc	pos|score	pos|p-value	pos|fdr	pos|rank	pos|goodsgrna	pos|lfc\n'
      'RPL19	4	2.695e-09	2.5851e-07	0.000707	3	4	-28.46707	1.0	1.0	1.0	19148	0	0\n'
      'KIF18B	5	1.0136e-08	2.5851e-07	0.000707	4	5	-26.55594	1.0	1.0	1.0	19146	0	0\n'
    )
    filepath_condition2 = '/folder/condition2.txt'
    self.fs.create_file(filepath_condition2, contents=file_content_condition2)

    files = [filepath_condition1, filepath_condition2]

    expected = pd.DataFrame([
      { 'condition': 'condition1', 'id': 'RPS5', 'num': 5, 'neg|score': 5.9353e-10, 'neg|p-value': 2.5851e-07, 'neg|fdr': 0.000707, 'neg|rank': 2, 'neg|goodsgrna': 5, 'neg|lfc': -30.64996, 'pos|score': 1.0, 'pos|p-value': 1.0, 'pos|fdr': 1.0, 'pos|rank': 19149, 'pos|goodsgrna': 0, 'pos|lfc': 0  },
      { 'condition': 'condition1', 'id': 'GTF2B', 'num': 5, 'neg|score': 2.0462e-10, 'neg|p-value': 2.5851e-07, 'neg|fdr': 0.000707, 'neg|rank': 1, 'neg|goodsgrna': 5, 'neg|lfc': -32.18633, 'pos|score': 1.0, 'pos|p-value': 1.0, 'pos|fdr': 1.0, 'pos|rank': 19150, 'pos|goodsgrna': 0, 'pos|lfc': 0  },
      { 'condition': 'condition2', 'id': 'RPL19', 'num': 4, 'neg|score': 2.695e-09, 'neg|p-value': 2.5851e-07, 'neg|fdr': 0.000707, 'neg|rank': 3, 'neg|goodsgrna': 4, 'neg|lfc': -28.46707, 'pos|score': 1.0, 'pos|p-value': 1.0, 'pos|fdr': 1.0, 'pos|rank': 19148, 'pos|goodsgrna': 0, 'pos|lfc': 0  },
      { 'condition': 'condition2', 'id': 'KIF18B', 'num': 5, 'neg|score': 1.0136e-08, 'neg|p-value': 2.5851e-07, 'neg|fdr': 0.000707, 'neg|rank': 4, 'neg|goodsgrna': 5, 'neg|lfc': -26.55594, 'pos|score': 1.0, 'pos|p-value': 1.0, 'pos|fdr': 1.0, 'pos|rank': 19146, 'pos|goodsgrna': 0, 'pos|lfc': 0  },
    ])
    self.assertEqual(convert_mageck(files), expected)

class ConvertRanks(pyfakefs.fake_filesystem_unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e

  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
    self.setUpPyfakefs()

  def test_v1(self):
    file_content_condition1 = (
      'Gene\tRANKS_score\tp-value\tFDR\t#_of_sgRNAs_considered\n'
      'EQUG\t-1.38\t0.015925\t0.0159\t8\n'
      'EQUG\t1.76\t0.193294\t0.9094\t8\n'
      'ISZY\t-1.26\t0.116286\t0.9975\t8\n'
    )
    filepath_condition1 = '/folder/condition1.txt'
    self.fs.create_file(filepath_condition1, contents=file_content_condition1)
    file_content_condition2 = (
      'Gene\tRANKS_score\tp-value\tFDR\t#_of_sgRNAs_considered\n'
      'VBEG\t3.77\t0.005852\t0.8344\t8\n'
      'VBEG\t-7.58\t0.1647\t0.0009\t8\n'
      'XWAS\t0.22\t0.168431\t0.0425\t8\n'
    )
    filepath_condition2 = '/folder/condition2.txt'
    self.fs.create_file(filepath_condition2, contents=file_content_condition2)

    files = [filepath_condition1, filepath_condition2]

    expected = pd.DataFrame([
      {
        'condition': 'condition1',
        'Gene': 'EQUG',
        'depletion_score': -1.38, 'depletion_p-value': 0.015925, 'depletion_FDR': 0.0159, 'depletion_#_of_sgRNAs_considered': 8,
        'enrichment_score': 1.76, 'enrichment_p-value': 0.193294, 'enrichment_FDR': 0.9094, 'enrichment_#_of_sgRNAs_considered': 8,
      },
      {
        'condition': 'condition1',
        'Gene': 'ISZY',
        'depletion_score': -1.26, 'depletion_p-value': 0.116286, 'depletion_FDR': 0.9975, 'depletion_#_of_sgRNAs_considered': 8,
        'enrichment_score': 0, 'enrichment_p-value': 1, 'enrichment_FDR': 1, 'enrichment_#_of_sgRNAs_considered': 0,
      },
      {
        'condition': 'condition2',
        'Gene': 'VBEG',
        'depletion_score': -7.58, 'depletion_p-value': 0.1647, 'depletion_FDR': 0.0009, 'depletion_#_of_sgRNAs_considered': 8,
        'enrichment_score': 3.77, 'enrichment_p-value': 0.005852, 'enrichment_FDR': 0.8344, 'enrichment_#_of_sgRNAs_considered': 8,
      },
      {
        'condition': 'condition2',
        'Gene': 'XWAS',
        'depletion_score': 0, 'depletion_p-value': 1, 'depletion_FDR': 1, 'depletion_#_of_sgRNAs_considered': 0,
        'enrichment_score': 0.22, 'enrichment_p-value': 0.168431, 'enrichment_FDR': 0.0425, 'enrichment_#_of_sgRNAs_considered': 8,
      },
    ])
    self.assertEqual(convert_ranks(files), expected)

  def test_v2(self):
    file_content_condition1 = (
      'Gene\tRANKS_score\tp-value\tFDR\t#_of_sgRNAs_considered\n'
      'EQUG\t-1.38\t0.015925\t0.0159\t8\n'
      'LTZA\t1.76\t0.193294\t0.9094\t8\n'
      'ISZY\t-1.26\t0.116286\t0.9975\t8\n'
    )
    filepath_condition1 = '/folder/condition1.txt'
    self.fs.create_file(filepath_condition1, contents=file_content_condition1)
    file_content_condition2 = (
      'Gene\tRANKS_score\tp-value\tFDR\t#_of_sgRNAs_considered\n'
      'VBEG\t3.77\t0.005852\t0.8344\t8\n'
      'KHGG\t-7.58\t0.1647\t0.0009\t8\n'
      'XWAS\t0.22\t0.168431\t0.0425\t8\n'
    )
    filepath_condition2 = '/folder/condition2.txt'
    self.fs.create_file(filepath_condition2, contents=file_content_condition2)

    files = [filepath_condition1, filepath_condition2]

    expected = pd.DataFrame([
      { 'condition': 'condition1', 'Gene': 'EQUG', 'RANKS_score': -1.38, 'p-value': 0.015925, 'FDR': 0.0159, '#_of_sgRNAs_considered': 8 },
      { 'condition': 'condition1', 'Gene': 'LTZA', 'RANKS_score': 1.76, 'p-value': 0.193294, 'FDR': 0.9094, '#_of_sgRNAs_considered': 8 },
      { 'condition': 'condition1', 'Gene': 'ISZY', 'RANKS_score': -1.26, 'p-value': 0.116286, 'FDR': 0.9975, '#_of_sgRNAs_considered': 8 },
      { 'condition': 'condition2', 'Gene': 'VBEG', 'RANKS_score': 3.77, 'p-value': 0.005852, 'FDR': 0.8344, '#_of_sgRNAs_considered': 8 },
      { 'condition': 'condition2', 'Gene': 'KHGG', 'RANKS_score': -7.58, 'p-value': 0.1647, 'FDR': 0.0009, '#_of_sgRNAs_considered': 8 },
      { 'condition': 'condition2', 'Gene': 'XWAS', 'RANKS_score': 0.22, 'p-value': 0.168431, 'FDR': 0.0425, '#_of_sgRNAs_considered': 8 },
    ])
    self.assertEqual(convert_ranks(files), expected)

class ExtractConditionFromFilename(unittest.TestCase):
  def test(self):
    file = 'folder/condition1.txt'
    expected = 'condition1'
    self.assertEqual(extract_condition_from_filename(file), expected)

class GetColumnNames(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()
  
  def test(self):
    file_contents = (
      'Gene\tBF\n'
      'GENEA\t10\n'
      'GENEB\t5\n'
      'GENEC\t1\n'
    )
    filepath = '/folder/condition1.txt'
    self.fs.create_file(filepath, contents=file_contents)

    expected = ['Gene', 'BF']
    self.assertEqual(get_column_names(filepath), expected)

class GetFiles(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()
  
  def test(self):
    self.fs.create_file('/folder/condition1.txt', contents='')
    self.fs.create_file('/folder/condition2.txt', contents='')
    self.fs.create_file('/folder/condition3.txt', contents='')

    folder = 'folder'

    expected = [
      'folder/condition1.txt',
      'folder/condition2.txt',
      'folder/condition3.txt',
    ]
    self.assertEqual(get_files(folder), expected)

class MergeAndAddConditions(pyfakefs.fake_filesystem_unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e

  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
    self.setUpPyfakefs()
  
  def test(self):
    file_content_condition1 = (
      'Gene\tBF\n'
      'GENEA\t10\n'
      'GENEB\t5\n'
      'GENEC\t1\n'
    )
    filepath_condition1 = '/folder/condition1.txt'
    self.fs.create_file(filepath_condition1, contents=file_content_condition1)
    file_content_condition2 = (
      'Gene\tBF\n'
      'GENEA\t8\n'
      'GENEB\t2\n'
      'GENED\t-3\n'
    )
    filepath_condition2 = '/folder/condition2.txt'
    self.fs.create_file(filepath_condition2, contents=file_content_condition2)

    files = [filepath_condition1, filepath_condition2]

    expected = pd.DataFrame([
      { 'condition': 'condition1', 'Gene': 'GENEA', 'BF': 10 },
      { 'condition': 'condition1', 'Gene': 'GENEB', 'BF': 5 },
      { 'condition': 'condition1', 'Gene': 'GENEC', 'BF': 1 },
      { 'condition': 'condition2', 'Gene': 'GENEA', 'BF': 8 },
      { 'condition': 'condition2', 'Gene': 'GENEB', 'BF': 2 },
      { 'condition': 'condition2', 'Gene': 'GENED', 'BF': -3 },
    ])
    self.assertEqual(merge_and_add_conditions(files), expected)

class MoveConditionColumn(unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e

  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
  
  def test(self):
    df = pd.DataFrame([
      { 'Gene': 'GENEA', 'BF': 10, 'condition': 'condition1' },
      { 'Gene': 'GENEB', 'BF': 5, 'condition': 'condition1' },
      { 'Gene': 'GENEC', 'BF': 1, 'condition': 'condition1' },
      { 'Gene': 'GENEA', 'BF': 8, 'condition': 'condition2' },
      { 'Gene': 'GENEB', 'BF': 2, 'condition': 'condition2' },
      { 'Gene': 'GENED', 'BF': -3, 'condition': 'condition2' },
    ])

    expected = pd.DataFrame([
      { 'condition': 'condition1', 'Gene': 'GENEA', 'BF': 10 },
      { 'condition': 'condition1', 'Gene': 'GENEB', 'BF': 5 },
      { 'condition': 'condition1', 'Gene': 'GENEC', 'BF': 1 },
      { 'condition': 'condition2', 'Gene': 'GENEA', 'BF': 8 },
      { 'condition': 'condition2', 'Gene': 'GENEB', 'BF': 2 },
      { 'condition': 'condition2', 'Gene': 'GENED', 'BF': -3 },
    ])
    self.assertEqual(move_condition_column(df), expected)
