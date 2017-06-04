#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test the search module
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import sys
import pytest
import re


@pytest.mark.parametrize(("pattern", 'mismatches', 'exphits_path'), [
    (b'TGGATGTGAAATGAGTCAAG', 3, 'tests/data/TGGATGTGAAATGAGTCAAG-results.sam'),
    (b'GGGTGGGGGGAGTTTGCTCC', 3, 'tests/data/vegfa-site1-results.sam'),
])
def test_search(pattern, mismatches, exphits_path):
    """Check if the search finds only (and all of) the expected hits."""
    # Need to explicitly add this path in case the user has installed 
    # big-search using the --user flag to setup install.
    sys.path.append(os.path.expanduser('~/.local/lib/python2.7/site-packages'))
    # Also import the module with my solution.
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import dictionarysearch as ds
    
    # If needed, preprocess the reference and save the resulting index.
    if (not os.path.isdir(ds.create_path(path_type='index_folder'))
        or not os.listdir(ds.create_path(path_type='index_folder'))):
        ds.preprocess_reference(ds.create_path(path_type='raw_ref'))
    
    # Get hits from my solution.  Note that SAM files (like the sample
    # files provided here) obtained using standard aligners report the
    # pattern sequence but not the reference sequence, so that will not
    # be parsed here.
    result = []
    with open(ds.search(pattern=pattern, max_mismatches=mismatches), 'rb') as myhits:
        for line in myhits:
            result.append(re.search(r'^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)', line).groups())
    
    # Parse file containing expected outputs.
    expected_hits = []
    with open(exphits_path, 'rb') as exphits:
        for line in exphits:
            # Unfortunately, pysam is not able to parse the 2 files with
            # the expected hits, as it complains they have corrupted
            # headers (sequence id lines are missing).  The easiest
            # workaround is to simply use re instead.
            if not line.startswith('@'):
                # Note: we also read the sequence, although that is not
                # necessary as it is determined uniquely by the second
                # field (flag field) in these simple cases.
                exphitRE = r'^\s*\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)'
                expected_hits.append(re.search(exphitRE, line).groups())
    
    # Check that my solution finds exactly the same hits.
    assert sorted(expected_hits) == sorted(result)
