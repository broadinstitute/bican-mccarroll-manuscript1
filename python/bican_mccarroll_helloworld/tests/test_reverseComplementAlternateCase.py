#!/usr/bin/env python3
# MIT License
# 
# Copyright 2025 Broad Institute
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import unittest
from bican_mccarroll_helloworld.reverseComplementAlternateCase import reverse_complement_alternate_case

class TestReverseComplementAlternateCase(unittest.TestCase):

    def test_reverse_complement_alternate_case(self):
        # Test cases
        test_cases = [
            ("ATCG", "CgAt"),
            ("atcg", "CgAt"),
            ("", ""),
            ("A", "T"),
            ("a", "T")
        ]

        for seq, expected in test_cases:
            with self.subTest(seq=seq):
                result = reverse_complement_alternate_case(seq)
                self.assertEqual(result, expected)