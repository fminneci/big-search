# big-search
Case Study in Python searching datasets larger than RAM.

Please find the detailed documentation for this project in the `docs` folder, 
specifically `docs/Report.md`.

A brief overview is outlined in the following instructions.

### Install

Clone this repository using:

    git clone https://github.com/fminneci/big-search.git

In the `big-search` directory, install with: 

    python setup.py install --user

...or with your favourite version of this command. As usual, this requires that 
you have `pip`/`setuptools` installed on you system.

Note that this solution uses Python's `pysam` module, as suggested. The module 
is installed automatically by the command above, however this requires zlib / 
libbzip2 / liblzma development files (they are used to install HTSlib within 
`pysam`). If needed, they can be obtained before running the above installation, 
for example with:
* `apt-get install zlib1g-dev libbz2-dev liblzma-dev`   (on Debian/Ubuntu Linux)
* `yum install zlib-devel bzip2-devel xz-devel`   (on RPM-based Linux or CygWin)

### Pre-processing and search steps

My solution is fully contained in the Python module named `dictionarysearch`.

This solution uses pre-processing of the reference genome (see the docs), along 
the lines of Example Strategy 3 - however, the reference index is based on an 
original idea rather than suffix arrays (this makes the search step simpler 
and faster, and use very little memory for the search).

Note that the reference genome needs to be then downloaded, using the link 
provided in the `ASSIGNMENT.md` file, to 
`tests/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`.

Running the module as a script with appropriate command line arguments will 
first pre-process the reference if needed, and then run the search.

This means that in order to pre-process the reference, you can simply run:

    python dictionarysearch.py TGGATGTGAAATGAGTCAAG 3

Note that this will also run the search at the end, which will save an output 
file named, in this case, `alignments_TGGATGTGAAATGAGTCAAG.txt` in the main 
directory.

For the pre-processing phase, this command will use 1 thread/process only, and 
roughly 500 MB memory, completing in 26 minutes on a machine with Intel i7 
(3.4 GHz). It will create a `ref_index` folder in the main directory, taking 
roughly 7 GB of disk space.

For the search phase, it will use multiprocessing (a maximum of 8 processes, and 
no more than the value returned by `multiprocessing.cpu_count()` on your machine). 
Depending on your system (number of processes, speed) the running time will vary 
considerably. On the same machine mentioned above, running 8 processes, the 
search will take around 68 seconds and use roughly 450 MB memory (around 60 MB 
memory per process). 

### Tests

After that, you can test that the search step of my solution passes all test as 
required (2 sample files provided).

Note that due to bad interaction with tox/pytest, one of their requirements (the 
`virtualenv` module) cannot be specified simply using the "tests_require" section 
of the `setup.py` file - it needs to be installed in advance. This can be achieved 
by running:

    pip install --user -r requirements.txt

...which also installs `pylint`.

Having done that, the tests can be performed by running:

    python setup.py test -a "-e py27"

This command will run tox/pytest (as usual, they will automatically install any 
missing necessary packages locally) and perform both tests, returning two 
passes.  

All tests will always include re-running the search step.

### Profiling

Finally, you can run the profiler on the search step, as required:

    python profilers/profile_search.py

This will use `timeit.repeat()` to run the search 3 consecutive times and report 
the running times in seconds - again, using roughly 60 MB memory per process.

