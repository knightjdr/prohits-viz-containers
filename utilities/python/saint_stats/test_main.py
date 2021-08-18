import pyfakefs.fake_filesystem_unittest
import unittest

from .main import read_saint, write_summary

class ReadSaint(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()
  
  def test(self):
    file_contents = (
      'Bait\tPrey\tPreyGene\t\t\t\t\t\t\t\t\t\t\t\tBFDR\n'
      'AAA_cond1-Nt\t111\tprey1\t\t\t\t\t\t\t\t\t\t\t\t0.0\n'
      'AAA_cond1-Nt\t222\tprey2\t\t\t\t\t\t\t\t\t\t\t\t0.01\n'
      'AAA_cond1-Nt\t444\tBBB\t\t\t\t\t\t\t\t\t\t\t\t0.01\n'
      'AAA_cond2-Ct\t111\tprey1\t\t\t\t\t\t\t\t\t\t\t\t0.0\n'
      'BBB_cond2-Ct\t111\tprey1\t\t\t\t\t\t\t\t\t\t\t\t0.01\n'
      'BBB_cond2-Ct\t222\tprey2\t\t\t\t\t\t\t\t\t\t\t\t0.01\n'
      'BBB_cond2-Ct\t333\tprey3\t\t\t\t\t\t\t\t\t\t\t\t0.02\n'
      'BBB_cond2-Ct\t555\tAAA\t\t\t\t\t\t\t\t\t\t\t\t0.01\n'
    )
    filepath = '/test/saint.txt'
    self.fs.create_file(filepath, contents=file_contents)

    fdr = 0.01

    expected = {
      'interactions': {
        'bait': {
          'AAA_cond1-Nt': 3,
          'AAA_cond2-Ct': 1,
          'BBB_cond2-Ct': 3,
        },
        'total': {
          'AAA_cond1-Nt-prey1': True,
          'AAA_cond1-Nt-prey2': True,
          'AAA_cond1-Nt-BBB': True,
          'AAA_cond2-Ct-prey1': True,
          'BBB_cond2-Ct-prey1': True,
          'BBB_cond2-Ct-prey2': True,
          'BBB_cond2-Ct-AAA': True,
        },
        'unique': {
          'AAA-prey1': True,
          'AAA-prey2': True,
          'AAA-BBB': True,
          'BBB-prey1': True,
          'BBB-prey2': True,
        },
      },
      'preys': {
        'significant': {
          'prey1': True,
          'prey2': True,
          'BBB': True,
          'AAA': True,
        },
        'total': {
          'prey1': True,
          'prey2': True,
          'prey3': True,
          'AAA': True,
          'BBB': True,
        },
      },
    }
    self.assertEqual(read_saint(filepath, fdr), expected)


class WriteSummary(pyfakefs.fake_filesystem_unittest.TestCase):
  def setUp(self):
    self.setUpPyfakefs()
  
  def test(self):
    self.fs.create_file('saint-statistics.txt', contents='')

    summary = {
      'interactions': {
        'bait': {
          'AAA_cond1-Nt': 3,
          'AAA_cond2-Ct': 1,
          'BBB_cond2-Ct': 3,
        },
        'total': {
          'AAA_cond1-Nt-prey1': True,
          'AAA_cond1-Nt-prey2': True,
          'AAA_cond1-Nt-BBB': True,
          'AAA_cond2-Ct-prey1': True,
          'BBB_cond2-Ct-prey1': True,
          'BBB_cond2-Ct-prey2': True,
          'BBB_cond2-Ct-AAA': True,
        },
        'unique': {
          'AAA-prey1': True,
          'AAA-prey2': True,
          'AAA-BBB': True,
          'BBB-prey1': True,
          'BBB-prey2': True,
        },
      },
      'preys': {
        'significant': {
          'prey1': True,
          'prey2': True,
          'BBB': True,
          'AAA': True,
        },
        'total': {
          'prey1': True,
          'prey2': True,
          'prey3': True,
          'AAA': True,
          'BBB': True,
        },
      },
    }

    write_summary(summary, 'saint.txt')

    expected = (
      'file: saint.txt\n'
      'interactions: unique - 5, total - 7\n'
      'preys: significant - 4, total - 5\n'
      'min: AAA_cond2-Ct with 1 prey(s)\n'
      'max: AAA_cond1-Nt with 3 prey(s)\n'
      'median: 3 prey(s)\n'
    )
    
    contents = ''
    with open('saint-statistics.txt') as f:
      contents = f.read()
    self.assertEqual(contents, expected)
