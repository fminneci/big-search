#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Naive exploration of search algo
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals


def compare_np5(large, small, mism):
    import numpy as np
    largenp = np.array([ord(c) for c in large], dtype='u1')
    smallnp = np.array([ord(c) for c in small], dtype='u1')
    S = smallnp.size
    out = []
    for i in xrange(largenp.size - S + 1):
        if np.count_nonzero(largenp[i:i+S] != smallnp) <= mism:
            out.append(i)
    return out


def compare_bin8(large, small, mism):
    mycode = {'A': b'00', 'C': b'01', 'G': b'10', 'T': b'11'}
    large_str = b''.join(mycode[c] for c in large)
    small_str = b''.join(mycode[c] for c in small)
    small_num = int(small_str, 2)
    S = len(small)
    out = []
    for i in range(len(large) - S + 1):
        large_num = int(large_str[2*i:2*(i+S)], 2)
        xor_num = small_num ^ large_num
        or_num = xor_num | 2 * xor_num
        if bin(or_num)[-2:1:-2].count('1') <= mism:
            out.append(i)
    return out


def compare_bin15(large, small, mism):
    import gmpy2
    mycode = {'A': b'1000', 'C': b'0100', 'G': b'0010', 'T': b'0001'}
    large_str = b''.join(mycode[c] for c in large)
    small_str = b''.join(mycode[c] for c in small)
    small_num = int(small_str, 2)
    S = len(small)
    out = []
    for i in range(len(large) - S + 1):
        large_num = int(large_str[4 * i : 4 * (i + S)], 2)
        and_num = small_num & large_num
        if gmpy2.popcount(and_num) >= S - mism:
            out.append(i)
    return out


def compare_np16(large, small, mism):
    import numpy as np
    L, S = len(large), len(small)
    largenp = np.array([ord(c) for c in large], dtype='u1')
    smallnp_init = np.empty(L, dtype='u1')
    for i in range(S):
        smallnp_init[i::S] = ord(small[i])
    rep1, extra = divmod(L, S)
    rep2 = rep1 - 1
    out = []
    for i in xrange(S):
        smallnp = np.roll(smallnp_init, i)
        xor_np = largenp != smallnp
        repetitions = rep1 if i <= extra else rep2
        for j in xrange(repetitions): 
            start = i + j * S
            if np.count_nonzero(xor_np[start:start+S]) <= mism:
                out.append(start)
    return sorted(out)


def seedextend_np5_re(large, small, mism):
    import re
    import numpy as np
    largenp = np.array([ord(c) for c in large], dtype='u1')
    smallnp = np.array([ord(c) for c in small], dtype='u1')
    L, S = len(large), len(small)
    r = S // (mism + 1)
    out = []
    for subsmall_offset in range(S - r + 1):
        subsmall = small[subsmall_offset:subsmall_offset+r]
        for submatch in re.finditer(subsmall, large[subsmall_offset:L-S+subsmall_offset+r]):
            start = submatch.start() # Because (due to line above) it is =   (submatch.start() + subsmall_offset) - subsmall_offset
            if start in out:
                continue
            if np.count_nonzero(largenp[start:start+S] != smallnp) <= mism:
                out.append(start)
    return sorted(out)


if __name__ == '__main__':
    
    mylong = b''.join(b'ACTG'[randint(0, 3)] for _ in range(10**6))
    myshort = b'TGGATGTGAAATGAGTCAAG'
    
    num_mismatches = 3
    
    for _ in range(10): 
        #print compare_np5(mylong, myshort, num_mismatches)
        #print compare_bin8(mylong, myshort, num_mismatches)
        #print compare_bin15(mylong, myshort, num_mismatches)
        print compare_np16(mylong, myshort, num_mismatches)
        #print seedextend_np5_re(mylong, myshort, num_mismatches)

