# Copyright (c) 2024 Joseph Hale
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import unittest

from run_batch import parse_args, parse_csv_line, ConfigFile


class Test_CLI(unittest.TestCase):

    def test_reads_options_correctly(self):
        parsed = parse_args(
            [
                "--csvpath",
                "batch.csv",
                "--rundir",
                "tycho_workspace",
                "--modeldir",
                "models",
            ]
        )
        self.assertEqual(parsed["csvpath"], "batch.csv")
        self.assertEqual(parsed["rundir"], "tycho_workspace")
        self.assertEqual(parsed["modeldir"], "models")

    def test_raises_when_csvpath_missing(self):
        with self.assertRaises(SystemExit):
            parse_args(["--rundir", "tycho_workspace", "--modeldir", "models"])

    def test_raises_when_rundir_missing(self):
        with self.assertRaises(SystemExit):
            parse_args(["--csvpath", "batch.csv", "--modeldir", "models"])

    def test_raises_when_modeldir_missing(self):
        with self.assertRaises(SystemExit):
            parse_args(["--csvpath", "batch.csv", "--rundir", "tycho_workspace"])


class Test_CSV_Parser(unittest.TestCase):

    def test_separates_columns_by_prefix(self):
        line = {
            "thing": "no-prefix-thing",
            "a.first": "a-first",
            "a.second": "a-second",
            "b.first": "b-first",
            "b.second": "b-second",
        }
        parsed = parse_csv_line(line)
        self.assertEqual(parsed[""], {"thing": "no-prefix-thing"})
        self.assertEqual(parsed["a"], {"first": "a-first", "second": "a-second"})
        self.assertEqual(parsed["b"], {"first": "b-first", "second": "b-second"})


GENEX_CONFIG = """
'.....................................................................'
' GENEX:                                                              '
'.....................................................................'
' Scaling(0=M,R;1=sol scale;2=OPAL)..........' 'ifm'       0
' lumin.ne.0..reevaluate luminosity=radiative' 'lumin'     0
' OPAL opacity type (0=type2,1=type1)........' 'nopac'     1
' OPAL EOS (0=no,1=yes)......................' 'nopaleos'  1
'..............................................MASS and RADIUS scaling'
' ifm=0:Fractional mass scaling   (0=ignore).' 'fm'        0.7d0
' ifm=0:Fractional radius scaling (0=ignore).' 'fr'        0.7d0
'.........................................SYSTEMATIC ABUNDANCE SCALING'
' ifm=1:metallicity scale relative to solar..' 'zscale'    0.0d-0
' ifm=1:production ratio He4/z...............' 'hetoz'     2.1d0
'......................................PARAMETERIZED ABUNDANCE SCALING'
' ifm=2:mass fraction for metallicity........' 'zpop'      0.01881d0
' ifm=2:mass fraction for He4................' 'xhe'       0.2676
' ifm=2:mass fraction for extra C12 (OPAL)...' 'xc12'      0.00
' ifm=2:mass fraction for extra O16 (OPAL)...' 'xo16'      0.00
'.....................................................................'
' Solid body rotation if omegbar .gt. 0/sec..' 'omegbar'   6.1d-8
'.....................................................................'
' Mixing mode flag...........................' 'mixmode'   2
'.....................................................................'
"""

class Test_ConfigFile(unittest.TestCase):

    def test_can_read_variables_like_a_dict(self):
        genex = ConfigFile(GENEX_CONFIG)
        self.assertEqual(genex["mixmode"], "2")
    
    def test_reading_a_nonexistent_variable_raises_KeyError(self):
        genex = ConfigFile(GENEX_CONFIG)
        with self.assertRaises(KeyError):
            genex["cookies"]

    def test_setting_a_nonexistent_variable_raises_KeyError(self):
        genex = ConfigFile(GENEX_CONFIG)
        with self.assertRaises(KeyError):
            genex["cookie"] = "monster"
        self.assertEqual(str(genex), GENEX_CONFIG)

    def test_setting_a_valid_variable_updates_the_config(self):
        genex = ConfigFile(GENEX_CONFIG)
        genex["xhe"] = "0.1234"
        self.assertEqual(str(genex), GENEX_CONFIG.replace("0.2676", "0.1234"))