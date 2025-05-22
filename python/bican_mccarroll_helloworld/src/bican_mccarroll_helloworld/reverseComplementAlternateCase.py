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
"""
A function that takes a string and returns its reverse complement with the case of the letters alternating
between upper and lower case.

This package depends on Biopython in order to demonstrate creating dependencies.
"""
from Bio.Seq import Seq

def reverse_complement_alternate_case(seq: str) -> str:
    """
    Returns the reverse complement of a DNA sequence with alternating case.

    :param seq: A string representing a DNA sequence.
    :return: A string representing the reverse complement of the input sequence with alternating case.
    """
    # Create a Seq object
    seq_obj = Seq(seq)

    # Get the reverse complement
    rev_comp = str(seq_obj.reverse_complement())

    # Alternate case
    alt_case = ''.join(
        char.upper() if i % 2 == 0 else char.lower() for i, char in enumerate(rev_comp)
    )

    return alt_case