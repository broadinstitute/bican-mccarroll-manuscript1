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
A simple hello world program that demonstrates how to use argparse to parse command line arguments,
and how to configure unit tests.
"""
import argparse
import sys

def hello(options):
    """
    A simple hello world function that returns a greeting message.
    :param options: The command line options.
    :return: A greeting message.
    """
    if options.mom:
        return "Hi, Mom!"
    else:
        return "Hello, World!"


def run(options):
    """
    Run the hello world program.
    :param options: The command line options.
    :return: 0 on success, 1 on failure.
    """
    print(hello(options))
    return 0

def main(args=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mom", "-m", default=False, action="store_true",
                        help="Say hello to mom instead of the world.")
    options = parser.parse_args(args)
    return run(options)

if __name__ == "__main__":
    sys.exit(main())
