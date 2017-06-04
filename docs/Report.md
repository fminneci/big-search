# Report on the big-search Python case study.

---

## My solution

The next three sections describe my solution to the big-search case study in 
detail.  
The first one describes the solution itself, contained in the `dictionarysearch` 
Python module; the following section describes tests and profiling steps taken; 
the third section summarises the conclusions.

### dictionarysearch

Considering that the task mentioned processing datasets larger than memory, and 
following Guido Van Rossum's linked post, it seemed important to obtain a 
solution that:
* searched using little memory
* explored employing some buffering of the input and output

Among the suggested strategies, the third one (saving a pre-processed reference 
index, based on suffix arrays or similar structures) was the most interesting 
to me - also see the Appendix for some early exploration of other routes. 
However, creating a pre-processed index does not help with buffering the input 
as the index structure is random access, thus it needs to be loaded all at once. 
Moreover, I wanted to implement multiprocessing in order to improve speed.

These considerations left an obvious route open: splitting the reference in 
small chunks during pre-processing, each with its own index, so that 
multiprocessing could act on these chunks with little memory usage during the 
search itself.  
Note that input buffering could still be used in the pre-processing step, but 
this seemed unnecessary as the pre-processing step is not the focus of the 
profiling/analysis.  
Buffering can still be used (and it is, in my solution) for the output.

Before starting a description of the module, one note about PEP8.  
I have written all of my code following Python idioms and PEP8 conventions. 
However, as there are two versions of the most controversial of those, I need to 
clarify that I have wrapped my code following the "72-99" version of the 
wrapping convention, rather than the "72-79" version. This means that my 
comments and docstrings are wrapped after 72 characters, while my code is 
wrapped after 99 characters. I thought this made for more understandable code.

As mentioned in the README file, my solution is fully contained in the Python 
module named `dictionarysearch`. It employs pre-processing of the reference 
genome, so that during the search only the pre-processed index (consisting of 
several files) is loaded, along the lines of Example Strategy 3 - however, the 
index is based on an original idea rather than suffix arrays.

The idea stems from the consideration that, because we need to account for 
mismatches, we cannot naively use the suffix array for the actual search - 
indeed, what we use at the "binary search" step in the outlined Example Strategy 
3 is only the initial character of every suffix, as we look up only one of the 
pattern's characters at a time. For this reason, we can avoid saving suffix 
arrays, which are complicated to prepare, and simply save a dictionary that 
associates each character in the reference's alphabet to a list of positions 
where that character occurs (basically, saving only the 1-character postings 
lists mentioned in Example Strategy 1). At the search stage, we won't even need 
to do the binary search: simply, for each position where the reference has the 
same character as the pattern character we are looking up (we get this information 
from the saved index dictionary) we increment the number of matches for that 
location. Finally, we look for location on the reference where the total number 
of matches is sufficient.

Obviously, there is no need to save an actual Python dictionary of lists - in 
the interest of saving memory, and along the lines of Example Strategy 3, we can 
simply save a number of numpy arrays in a single compressed numpy archive for 
each chunk of the reference. These will be loaded at the search stage in a way 
that is easy to parallelise using the `multiprocessing` module.

This can be seen in the code within the `dictionarysearch` module.

After the imports, the chunk size and the overlap between consecutive chunks 
(this is needed to avoid missing hits halfway between chunks) are defined on 
lines 20-22. Note that CHUNK_OVERLAP is also the shortest possible size for a 
chunk, thanks to line 57.

Next, an auxiliary function is defined (lines 25-40) to have consistent naming 
of files and directories.

On lines 43-77, the function for pre-processing the reference is defined - this 
will use pysam to open one contig (chromosome or short contig) at a time, split 
it in chunks (l. 54-59), save the chunk's sequence (l. 61-66), create and save 
the above mentioned numpy arrays in compressed format (l. 69-77).

Next is a helper function that returns the DNA reverse complement of a pattern 
sequence (as the reverse strand of the genome will need to be searched, too).

Lines 85-156 perform the search. Note that 3 out of the 4 functions defined here 
are purely needed to implement multiprocessing using Python's 
`multiprocessing.Pool` object. The actual search is done by `search_chunk`, 
which loads the reference sequence (l. 93-95) and the index arrays into a 
dictionary (l. 97-99), counts the number of matches for a pattern aligned at 
each position in the reference (l. 104-106), extract the starting positions with 
appropriate number of matches (l. 111-112), and returns a 5-tuple for each match 
(117-122).

I'd like to acknowledge the fact that the idea of using a numpy array of keeping 
track of the number of matches/mismatches, rather than a list of tuples or slower 
structures, was already used in another solution by user "chingtoe365". However, 
incrementing only the matching positions, rather than the mis-matching positions, 
improves the speed significantly.

Note that the `match_count` array is defined (l. 104) with data type appropriate 
for 1-byte integers: this means that the pattern cannot be longer than 256 
characters, which is anyway the case for short genomic "reads".

Each hit in the returned list of 5-tuples contains:

0. binary flag, either 0 (meaning "found on reference forward strand") or 16 
("on reverse strand")
1. name of the reference contig/chromosome
2. 1-based position of start of match on the reference
3. pattern sequence (either original or reverse complement)
4. corresponding sequence on the reference

This is because the output will be compared to the sample output files (standard 
genomics alignment files in SAM format), which have several defined features:
* positions are 1-based
* both the pattern sequence, and a flag to indicate forward/reverse strand, are 
reported

Moreover, the reference sequence is reported in my solution because this is how 
I interpreted one of the Guidelines (number 6) in the `ASSIGNMENT.md` file.

Finally, note that the Pool instance's `imap_unordered` method would itself write 
solutions to file as it receives them from the workers, therefore proving 
automatically a sort of buffer mechanisms - however, I implemented a proper 
buffer on top of that (l. 145-155).

The remaining lines cover the `main()` function, performing a few checks on the 
input parameters.

### Benchmarking

All benchmarking has been done using two tools:
* Python's `timeit` module, for monitoring running time
* saving output of Unix `top` utility, for monitoring memory usage

The benchmarking was performed on a machine with HDD disk, DDR3 memory, Intel i7 
quad-core (max 8 threads) with processor speed 3.4 GHz, running a Linux OS 
(Ubuntu 16.04).

For the test suite, it is easy to run it as indicated in the README file, 
obtaining a full pass outcome (both sample cases).

For the profiling script, again see the README file - only the running time is 
monitored with the profiling script, however I have separately recorded memory 
usage and timing for the various sub-sections of the code (data not shown).

As required, the profiling includes reading pattern and max number of mismatches 
from command line, loading the reference index, performing the search, and 
writing results to file.

One important caveat is that when running the search with multiprocessing the 
first time, processes were not used optimally and the search time was higher 
than later searches. I have not been able to investigate this properly, but the 
benchmarking is done on the later searches, as the initial slower run was badly 
reproducible (it slowed down the search by a variable factor, and the culprit 
seemed to be slower index loading times).

The detailed profiling data mentioned above indicate that:
* the pre-processing step (1 process, 500 MB memory, 26 minutes, 7 GB output 
index folder) had two time consuming steps, creating the arrays (27% time) and 
serialising/compressing/saving them to disk (70% time), with all other steps 
(initialising the arrays and counting the occurrences of different characters, 
loading, splitting, and saving the sequence files) adding up to 3% of the time
* the search step (8 processes, 60 GB memory per process, 68 seconds) had two 
time consuming steps, loading the arrays and sequences (31% time) and searching 
(63% time), while all other steps (initialising variables and local arrays, 
collecting results, printing them to file) added up to 6% of the time - in 
particular, the difference when varying the output buffer size on line 151 was 
negligible due to very small size of the output list

The sizable delay due to saving or loading the arrays might be reduced 
using a faster disk type (SSD), but I could not test that.

Note that the size of pre-processed reference chunks can be varied. After 
testing several sizes, I settled for one that gave low memory usage on 8 
processes as well as number of files in the reference index (npz files) smaller 
than 1000.

I also carried out further exploration, as follows.

Varying the values of _k_ and len(_P_), no change was detected in memory usage, 
while time varied as follows:

     k   len(P) : time (s)
    -----------------------
     1     16   :   57
     2     16   :   57
     3     16   :   57
    
     2     18   :   63
     3     18   :   63
     4     18   :   63
    
     2     20   :   68
     3     20   :   68
     4     20   :   68
    
     6     30   :   98
     7     30   :   98
    10     30   :   99
    
    15     50   :  159
    20     50   :  159

There is no real dependence on _k_, while times increase linearly with len(_P_).

This is expected, as the search step (lines 105-106) is clearly of time complexity 
O(len(_P_)), while it does not depend on _k_, which is only a threshold.

### Conclusions

I really enjoyed working on this project, as I found it is positioned exactly at 
the ideal intersection for me between programming/computing work (most of it) 
and bio-medical application (smaller part at the end, but significant and 
essential to complete the task fully), and I have been able to exploit all of my 
background.

In summary, I found that my solution has the following strong points:
* it is completely general on _k_, len(_P_), len(alphabet)
* running time is completely independent of k
* it was the fastest in my tests

However, running time is much slower than production levels of less than 1 
second. This is hopefully compensated by a very compact solution (one module with 
everything in it, 177 lines of code including many for comments or auxiliary 
functions).

A weak point perhaps would be the fact that running time scales with O(len(_P_)). 
However, for the length under examination in this task (around 20 characters), 
this is still the fastest solution I found.  
If used for much longer patterns, possibly a different solution can be used (for 
example a fast on-line seed-and-extend search like the one mentioned in the 
Appendix).

Also, possible improvements could be:
* depending on the task, the reference sequence for the chunks may not be not 
necessary - for example, it is not necessary if the output needs to be 
comparable to standard aligners' SAM files - and in this case these ".seq" files 
do not need to be created at all: the index arrays alone will do the job at the 
search stage
* the many index files could be stored in a more compact way
* command line parameters can become STDIN inputs if needed
* the number of used processes or the output buffer could be made into further 
parameters provided by the user

---

## Appendix

This Appendix outlines some of the preliminary exploration of the problem that 
I carried out before getting to my solution.

These can be found in an unpolished form in the accompanying script named 
`exploration_search.py`, which contains a few algorithms I wrote and tested 
initially.  
They are simply in the form of top-level functions, accepting a long string 
("large", i.e. the text), a short string ("small", i.e. the pattern), and a max 
number of mismatches as input parameters, performing the search, and returning 
the positions only. In this unrefined version, they do not implement 
multiprocessing nor any kind of sophisticated sequence input/output: for these 
preliminary tests, they were simply fed a large text string from memory.

These algorithms use either numpy arrays of integer-transformed version of the 
text, or binary encodings of the text, to do naive on-line searching. I was 
testing these in alternation using Unix `time` on the command line.

Note that `compare_bin15` requires install the module `gmpy2` module, which is 
not included here. In fact, these algorithms are not meant to be used - they are 
shown here only for reference purpose.

Algorithms `compare_np5`, `compare_bin8`, `compare_bin15`, `compare_np16` were 
the fastest naive on-line searching algorithms I obtained (`compare_np16` was 
the fastest); `seedextend_np5_re` was the seed-and-extend version of one of 
those, being roughly 3X times faster than the others.  
All of these were however slower than my solution described above, for patterns 
similar to the ones in the sample cases - note that the functions described here 
do not even have to load the sequence from disk.

