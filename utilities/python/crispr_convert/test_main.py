import pandas as pd
import pandas.testing as pd_testing
import pyfakefs.fake_filesystem_unittest
import unittest

from .main import (
  convert_bagel,
  extract_condition_from_filename,
  get_files,
  move_condition_column,
)

class ExtractConditionFromFilename(unittest.TestCase):
  def test(self):
    file = 'folder/condition1.txt'
    expected = 'condition1'
    self.assertEqual(extract_condition_from_filename(file), expected)

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

class ConvertBagel(pyfakefs.fake_filesystem_unittest.TestCase):
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
    self.assertEqual(convert_bagel(files), expected)
