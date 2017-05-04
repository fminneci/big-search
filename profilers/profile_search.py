#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Profile the main search module
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import sys
import timeit


def profile_genome_search():
    """Profile the search module."""
    #genome_file  = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
    #pattern = b'TGGATGTGAAATGAGTCAAG'
    #maxmismatch = 3
    
    # If needed, preprocess the reference file and save the resulting
    # index.  This section is not included in the profiling, as per the
    # instructions.
    if (not os.path.isdir(ds.create_path(path_type='index_folder'))
        or not os.listdir(ds.create_path(path_type='index_folder'))):
        ds.preprocess_reference(ds.create_path(path_type='raw_ref'))
    
    # Show on STDOUT three measurements of the duration of the full
    # search in seconds, including finding the alignments and printing
    # results to disk.
    statement = """ds.search(pattern=pattern, max_mismatches=maxmismatch)"""
    setup = """import dictionarysearch as ds; pattern = b'TGGATGTGAAATGAGTCAAG'; maxmismatch = 3"""
    print timeit.repeat(stmt=statement, setup=setup, number=1)


def main():
    profile_genome_search()


if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    import dictionarysearch as ds
    main()
