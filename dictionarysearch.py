#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Multiprocessing dictionary-like solution of the big-search problem.
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import sys
import glob
import re
from multiprocessing import cpu_count, Pool
from string import maketrans

import pysam
import numpy as np

SEQFILE_EXT = 'seq'  # Extension for processed reference sequence files
CHUNK_SIZE = 5 * 10**6  # Desired max size for a contig part, or "chunk"
CHUNK_OVERLAP = 1000  # Sequence overlap of consecutive chunks


def create_path(path_type, identifier=None):
    """Helper function for creating file paths."""
    main_dir = os.path.dirname(os.path.abspath(__file__))
    reference_dir = os.path.join(main_dir, 'ref_index')
    if path_type == 'raw_ref':
        return os.path.join(main_dir, 'tests/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
    elif path_type == 'proc_ref':
        return os.path.join(reference_dir, 'GRCh38_{0}.{1}'.format(identifier, SEQFILE_EXT))
    elif path_type == 'index':
        return os.path.join(reference_dir, 'GRCh38_{0}.index.npz'.format(identifier))
    elif path_type == 'index_folder':
        return reference_dir
    elif path_type == 'alignments':
        return os.path.join(main_dir, 'alignments_{0}.txt'.format(identifier))
    else:
        raise ValueError('Path type was not recognised - aborting!')


def preprocess_reference(ref_file):
    """Pre-process the genome reference."""
    # Make sure reference genome has been downloaded.
    if not os.path.exists(ref_file):
        raise IOError('Reference file not found - please download the reference genome')
    # If needed, create a folder for the reference index files.
    if not os.path.isdir(create_path(path_type='index_folder')):
        os.mkdir(create_path(path_type='index_folder'))
    # Loop through the reference contigs.
    with pysam.FastxFile(ref_file) as raw_ref:
        for contig in raw_ref:
            chunk = covered = start = 0
            end = len(contig.sequence)
            # Create overlapping chunks - no smaller than CHUNK_OVERLAP.
            while covered < end:
                chunkname = 'CONTIG_{0}_CHUNK_{1}'.format(contig.name, chunk)
                chunkseq = contig.sequence[start : start + CHUNK_SIZE + CHUNK_OVERLAP]
                # Save sequence chunks.
                with open(create_path(path_type='proc_ref', identifier=chunkname), 'wb') as ref:
                    ref.write(chunkseq)
                chunk += 1
                covered += CHUNK_SIZE + CHUNK_OVERLAP
                start += CHUNK_SIZE
                contig_counts = {c: chunkseq.count(c) for c in set(chunkseq)}
                # Initialise a dictionary of arrays and fill them with
                # the positions for each base in hte reference chunk.
                (arrays, positions) = ({}, {})
                for c, count in contig_counts.iteritems():
                    arrays[c] = np.empty(count, dtype='i4')
                    positions[c] = 0
                for p, c in enumerate(chunkseq):
                    arrays[c][positions[c]] = p
                    positions[c] += 1
                # Save index array in compressed format.
                np.savez_compressed(create_path(path_type='index', identifier=chunkname), **arrays)


def DNA_reverse_complement(p):
    """Return the DNA reverse complement of a sequence."""
    return p[::-1].translate(maketrans(b'ACGT', b'TGCA'))


def search_chunk(pattern, max_mismatches, chunkname):
    """Search for pattern in a specific chunk of the reference index."""
    len_pattern = len(pattern)
    revpattern = DNA_reverse_complement(pattern)
    chunk = int(re.search(r'_CHUNK_(\d+)$', chunkname).group(1))
    contig_name = re.search(r'^CONTIG_(\S+)_CHUNK_', chunkname).group(1)
    positions_offset = CHUNK_SIZE * chunk
    # Load reference chunk sequence.
    with open(create_path(path_type='proc_ref', identifier=chunkname), 'rb') as ref:
        ref_sequence = ref.read()
    len_ref = len(ref_sequence)
    # Load all needed index arrays into a dictionary.
    with np.load(create_path(path_type='index', identifier=chunkname)) as index:
        arrays = {c: index[c] if c in index else np.empty(0, dtype='i4')
                  for c in set(pattern).union(set(revpattern))}
    chunk_hits = []
    for i, query in enumerate((pattern, revpattern)):
        # Count positive matches for all substrings of ref_sequence,
        # storing the number of matches at the first letter's position.
        match_count = np.zeros(len_ref, dtype='i1')
        for p, c in enumerate(query):
            match_count[arrays[c] - p] += 1
        # Find positions with a sufficient number of matches - excluding
        # the last (len_pattern - 1) positions, where alignment is not
        # possible.  Note that anyway, such positions were corrupted by
        # the assignment on the previous line, whenever p > 0.
        positions_array = np.nonzero(match_count[: len_ref-len_pattern+1]
                                     >= len_pattern - max_mismatches)[0]
        # Add hits to the output list.  Each hit is a 5-tuple - see the
        # accompanying docs for details.  Also, note that, for chunks
        # after the first one, the initial positions have been searched
        # already in the previos chunk.
        start = (CHUNK_OVERLAP - len_pattern + 1) if chunk else 0
        for pos in positions_array:
            if pos >= start:
                chunk_hits.append(((0, 16)[i], contig_name, positions_offset + pos + 1,
                                   query, ref_sequence[pos : pos+len_pattern]))
    return chunk_hits


def search_chunk_worker(chunkname):
    """Single-parameter version of search_chunk, for pool.map()."""
    return search_chunk(search_chunk_worker.pattern, search_chunk_worker.max_mismatches, chunkname)


def init_pool_worker(pattern, max_mismatches):
    """Helper function to initialise a pool worker."""
    search_chunk_worker.pattern = pattern
    search_chunk_worker.max_mismatches = max_mismatches


def search(pattern, max_mismatches):
    """Search for pattern using the previously saved reference index."""
    filepaths = glob.iglob(os.path.join(create_path(path_type='index_folder'), '*' + SEQFILE_EXT))
    chunknames = (re.search(r'(CONTIG_\S+_CHUNK_\d+)', path).group(1) for path in filepaths)
    n_proc = min(8, cpu_count())
    pool = Pool(processes=n_proc, initializer=init_pool_worker, initargs=(pattern, max_mismatches))
    # Distribute chunks to processes, then write outputs to file (for
    # the sake of simplicity, not quite a full SAM file) as they are
    # received from the workers, but with a buffer.
    outpath = create_path(path_type='alignments', identifier=pattern)
    outlist = []
    print_string = r'{0:>3} {1:>12} {2:>12} {3:>24} {4:>24}'
    with open(outpath, 'wb', 0) as outfile:
        for result in pool.imap_unordered(search_chunk_worker, chunknames, 5):
            outlist.extend(result)
            if len(outlist) >= 20:
                outfile.write('\n'.join(print_string.format(*hit) for hit in outlist) + '\n')
                outlist = []
        if outlist:
            outfile.write('\n'.join(print_string.format(*hit) for hit in outlist) + '\n')
    return outpath


def main():
    # Get command-line parameters and check them.
    if len(sys.argv) != 3:
        raise TypeError('This program takes 2 parameters: pattern and max number of mismatches')
    _pattern, _max_mismatches = bytes(sys.argv[1]), int(sys.argv[2])
    if len(_pattern) > CHUNK_OVERLAP:
        raise ValueError('Pattern cannot be longer than {0} bp'.format(CHUNK_OVERLAP))
    
    # Preprocess the reference file, if needed.
    if (not os.path.isdir(create_path(path_type='index_folder'))
        or not os.listdir(create_path(path_type='index_folder'))):
        preprocess_reference(create_path(path_type='raw_ref'))
    # Find all alignments and print them to file.
    results_file = search(pattern=_pattern, max_mismatches=_max_mismatches)
    print 'Results successfully written to file: {0}'.format(results_file)


if __name__ == '__main__':
    main()
