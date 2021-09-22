import os
import pandas as pd
import pandas.testing as pd_testing
import pyfakefs.fake_filesystem_unittest
import unittest

from os import access

from .main import (
  check_file_in_excel_format,
  fix_data,
  read_file,
  read_symbols_to_fix,
  write_fixed_file,
  write_summary,
)

class CheckFileInExcelFormat(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()

  def test_text_file(self):
    file_contents = (
      'column1\tcolumn2\tcolumn3\n'
      '1-Sep\t111\tAAA\n'
      '1-Sep\t111\tBBB\n'
      '1-Sep\t111\tCCC\n'
    )
    filepath = '/test/file.txt'
    self.fs.create_file(filepath, contents=file_contents)

    self.assertFalse(check_file_in_excel_format(filepath))

  def test_xlsx_file(self):
    df = pd.DataFrame([
      { 'Bait': '1-Sep', 'Prey': '111', 'PreyGene': 'AAA' },
      { 'Bait': '1-Sep', 'Prey': '222', 'PreyGene': 'BBB' },
      { 'Bait': '1-Sep', 'Prey': '333', 'PreyGene': 'CCC' },
      { 'Bait': '1-Sep', 'Prey': '444', 'PreyGene': 'DDD' },
      { 'Bait': '1-Sep', 'Prey': '555', 'PreyGene': 'EEE' },
      { 'Bait': '1-Sep', 'Prey': '666', 'PreyGene': 'FFF' },
      { 'Bait': 'MARCH5', 'Prey': '111', 'PreyGene': 'AAA' },
      { 'Bait': 'MARCH5', 'Prey': '222', 'PreyGene': 'BBB' },
      { 'Bait': 'MARCH3', 'Prey': '777', 'PreyGene': 'GGG' },
    ])
    filepath = '/file.xlsx'
    df.to_excel(filepath)

    self.assertTrue(check_file_in_excel_format(filepath))

class ReadFile(pyfakefs.fake_filesystem_unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e
  
  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
    self.setUpPyfakefs()

  def test_text_file(self):
    is_excel = False
    file_contents = (
      'column1\tcolumn2\tcolumn3\n'
      '1-Sep\t111\tAAA\n'
      '1-Sep\t222\tBBB\n'
      '1-Sep\t333\tCCC\n'
    )
    filepath = '/test/file.txt'
    self.fs.create_file(filepath, contents=file_contents)

    expected = {
      'sheet': pd.DataFrame([
        { 'column1': '1-Sep', 'column2': 111, 'column3': 'AAA' },
        { 'column1': '1-Sep', 'column2': 222, 'column3': 'BBB' },
        { 'column1': '1-Sep', 'column2': 333, 'column3': 'CCC' },
      ])
    }

    actual = read_file(filepath, is_excel)
    self.assertEqual(actual['sheet'], expected['sheet'])

class ReadSymbolsToFix(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()

  def test(self):
    file_contents = (
      'Excel representation\tsymbol\tconverted\tofficial symbol\tambiguous\n'
      '5-Mar\tMARCH5\tMARCHF5\tFalse\n'
      '1-Sep\tSEP1\tXRN1\tTrue\n'
      '1-Sep\tSEPT1\tSEPTIN1\tTrue\n'
    )
    dir = os.path.abspath(os.path.dirname(__file__))
    filepath = f'{dir}/symbols-to-fix.txt'
    self.fs.create_file(filepath, contents=file_contents)

    expected = {
      '5-Mar': {
        'ambiguous': False,
        'converted_symbols': ['MARCH5'],
        'official_symbols': ['MARCHF5'],
      },
      '1-Sep': {
        'ambiguous': True,
        'converted_symbols': ['SEP1', 'SEPT1'],
        'official_symbols': ['SEPTIN1', 'XRN1'],
      },
    }

    self.assertEqual(read_symbols_to_fix(), expected)

class FixData(pyfakefs.fake_filesystem_unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e
  
  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

  def _test_fixed_data(self, arg0):
    columns = arg0
    data = {
      'sheet': pd.DataFrame([
        { 'column1': '2001-09-01', 'column2': '111', 'column3': 'AAA' },
        { 'column1': '1-Sep', 'column2': '222', 'column3': 'BBB' },
        { 'column1': '1-Sep', 'column2': '333', 'column3': '5-Mar' },
        { 'column1': '4-Oct', 'column2': '333', 'column3': '2021-03-05' },
        { 'column1': '4-Oct', 'column2': '333', 'column3': '2021-03-05 00:00:00' },
      ])
    }
    symbols_to_fix = {
      '5-Mar': {
        'ambiguous': False,
        'converted_symbols': ['MARCH5'],
        'official_symbols': ['MARCHF5'],
      },
      '4-Oct': {
        'ambiguous': False,
        'converted_symbols': ['Oct4'],
        'official_symbols': ['POU5F1'],
      },
      '1-Sep': {
        'ambiguous': True,
        'converted_symbols': ['SEP1', 'SEPT1'],
        'official_symbols': ['XRN1', 'SEPTIN1'],
      },
    }

    expected_fixed = {
      'sheet': pd.DataFrame([
        { 'column1': '1-Sep', 'column2': '111', 'column3': 'AAA' },
        { 'column1': '1-Sep', 'column2': '222', 'column3': 'BBB' },
        { 'column1': '1-Sep', 'column2': '333', 'column3': 'MARCHF5' },
        { 'column1': 'POU5F1', 'column2': '333', 'column3': 'MARCHF5' },
        { 'column1': 'POU5F1', 'column2': '333', 'column3': 'MARCHF5' },
      ])
    }
    expected_summary = [
      ('sheet', 'column1', '1-Sep'),
      ('sheet', 'column1', '4-Oct'),
      ('sheet', 'column3', '5-Mar'),
    ]
    actual_fixed, actual_summary = fix_data(data, columns, symbols_to_fix)
    self.assertEqual(actual_fixed['sheet'], expected_fixed['sheet'])
    self.assertEqual(actual_summary, expected_summary)

  def test_columns(self):
    self._test_fixed_data('column1|column3')

  def test_no_columns(self):
    self._test_fixed_data('')

class WriteSummary(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()

  def _test_written_summary(self, is_excel, expected):
    summary = [
      ('sheet1', 'column1', '1-Sep'),
      ('sheet1', 'column3', '5-Mar'),
    ]
    symbols_to_fix = {
      '5-Mar': {
        'ambiguous': False,
        'converted_symbols': ['MARCH5'],
        'official_symbols': ['MARCHF5'],
      },
      '1-Sep': {
        'ambiguous': True,
        'converted_symbols': ['SEP1', 'SEPT1'],
        'official_symbols': ['XRN1', 'SEPTIN1'],
      },
    }

    write_summary(summary, is_excel, symbols_to_fix)
    with open('./summary.txt', 'r') as f:
      actual = f.read()
    self.assertEqual(actual, expected)

  def test_is_excel(self):
    expected = (
      'sheet\tcolumn\toriginal symbol\tfixed\tcorrect(ed) symbol(s)\n'
      'sheet1\tcolumn1\t1-Sep\tFalse\tXRN1, SEPTIN1\n'
      'sheet1\tcolumn3\t5-Mar\tTrue\tMARCHF5\n'
    )

    self._test_written_summary(True, expected)

  def test_is_not_excel(self):
    expected = (
      'column\toriginal symbol\tfixed\tcorrect(ed) symbol(s)\n'
      'column1\t1-Sep\tFalse\tXRN1, SEPTIN1\n'
      'column3\t5-Mar\tTrue\tMARCHF5\n'
    )

    self._test_written_summary(False, expected)

class WriteFixedFile(pyfakefs.fake_filesystem_unittest.TestCase):
  def assertDataframeEqual(self, a, b, msg):
    try:
      pd_testing.assert_frame_equal(a, b)
    except AssertionError as e:
      raise self.failureException(msg) from e
  
  def setUp(self):
    self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
    self.setUpPyfakefs()

  def test_is_excel(self):
    data = {
      'sheet1': pd.DataFrame([
        { 'column1': '1-Sep', 'column2': 111, 'column3': 'AAA' },
        { 'column1': '1-Sep', 'column2': 222, 'column3': 'BBB' },
        { 'column1': '1-Sep', 'column2': 333, 'column3': 'MARCHF5' },
      ]),
      'sheet2': pd.DataFrame([
        { 'column1': 'POU5F1', 'column2': 111, 'column3': 'AAA' },
      ]),
    }
    filepath = 'file.xlsx'
    is_excel = True

    expected = {
      'sheet1': pd.DataFrame([
        { 'column1': '1-Sep', 'column2': 111, 'column3': 'AAA' },
        { 'column1': '1-Sep', 'column2': 222, 'column3': 'BBB' },
        { 'column1': '1-Sep', 'column2': 333, 'column3': 'MARCHF5' },
      ]),
      'sheet2': pd.DataFrame([
        { 'column1': 'POU5F1', 'column2': 111, 'column3': 'AAA' },
      ]),
    }

    write_fixed_file(data, filepath, is_excel)
    actual = pd.read_excel('./file-fixed.xlsx', sheet_name=None)
    self.assertEqual(actual['sheet1'], expected['sheet1'])
    self.assertEqual(actual['sheet2'], expected['sheet2'])

  def test_is_not_excel(self):
    data = {
      'sheet': pd.DataFrame([
        { 'column1': '1-Sep', 'column2': '111', 'column3': 'AAA' },
        { 'column1': '1-Sep', 'column2': '222', 'column3': 'BBB' },
        { 'column1': '1-Sep', 'column2': '333', 'column3': 'MARCHF5' },
      ])
    }
    filepath = 'file.txt'
    is_excel = False

    expected = (
      'column1\tcolumn2\tcolumn3\n'
      '1-Sep\t111\tAAA\n'
      '1-Sep\t222\tBBB\n'
      '1-Sep\t333\tMARCHF5\n'
    )

    write_fixed_file(data, filepath, is_excel)
    with open('./file-fixed.txt', 'r') as f:
      actual = f.read()
    self.assertEqual(actual, expected)