Protein back translation with pairwise distance constraints
============================================================

Usage
``````
Run using `python3 -m src.backTrans`

::

	usage: backTrans.py [-h] (-i INPUT_FILE | -s SEQUENCE) -m {dna,protein}
                    [-d MIN_DIST] [-n NUM_SAMPLES] [-o OUTPUT]
                    [--output_prefix OUT_PREFIX] [--stats_header]
                    [--forbidden FORBIDDEN]

    Backtranslate amino acid sequences with pairwise distance constraints

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT_FILE, --input_file INPUT_FILE
                            Path to fasta file containing amino acid or DNA
                            sequences
      -s SEQUENCE, --sequence SEQUENCE
                            amino acid sequence passed on command-line
      -m {dna,protein}, --mode {dna,protein}
                            mode to run the tool on: in protein mode, samples
                            backtranslations, in dna mode, uses the dna sequences
                            as samples directly.
      -d MIN_DIST, --min_distance MIN_DIST
                            Minimum distance (Hamming, as fraction of distinct
                            nucleotides) between all returned sequences
      -n NUM_SAMPLES, --num_samples NUM_SAMPLES
                            Maximum number of backtranslated DNA samples to
                            produce
      -o OUTPUT, --output-dir OUTPUT
                            An existing directory where the output will go. This
                            is a fasta of the sampled sequence(s) and a stats file
      --output_prefix OUT_PREFIX
                            prefix for output files
      --stats_header        Prints the header for the stats file
      --forbidden FORBIDDEN
                            File path to DNA sequences that cannot appear in the
                            sample; one sequence per line.


Protein(s)
-----------
Produces DNA sequences compatible with each protein, where no two DNA sequences have less than `--min_distance` Hamming distance.

* Provide a single amino acid sequence on command-line: `-s`
* Or an fasta file with one or more protein sequence: `-i` and `-m protein`

DNA(s)
-------
Finds a set of DNA sequences where no two have less than ``min--distance`` Hamming distance

Provide a fasta file: `-i` and `-m dna`


TODOs
``````

* Try CD-Hit on DNA sequences as opposed to graph approach. Note:

	* CD-Hit produces one representative per cluster, and the pairwise distance is only guaranteed lower than threshold on those
	* cd-hit-est seems to be capped at 75% ID minimum

* Faster pairwise distance computation:
	* cd-hit similar approach: avoid alignment for sequences where given number of shared small word frequency counts is exceeded.
	* sketching-based approach, eg Mash. Con: inexact
