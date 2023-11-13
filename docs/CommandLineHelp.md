# Command-Line Help for `sketchm`

This document contains the help content for the `SketchM` command-line program.

**Command Overview:**

* [`sketchm`↴](#sketchm)
* [`sketchm sketch`↴](#sketchm-sketch)
* [`sketchm index`↴](#sketchm-index)
* [`sketchm dist`↴](#sketchm-dist)
* [`sketchm info`↴](#sketchm-info)

## `sketchm`

Fast genome ANI estimation using k-mer sketches.

**Usage:** `sketchm [COMMAND]`

###### **Subcommands:**

* `sketch` — Create k-mer sketches for genomes
* `index` — Create k-mer index from sketches
* `dist` — Compute distances between genome sketches
* `info` — Display information about sketch file


## `sketchm sketch`

Create k-mer sketches for genomes

**Usage:** `sketchm sketch [OPTIONS] --output-file <OUTPUT_FILE> <--genome-files <GENOME_FILES>...|--genome-path-file <GENOME_PATH_FILE>>`

###### **Options:**

* `-f`, `--genome-files <GENOME_FILES>` — Genome FASTA/Q file(s) to sketch
* `-p`, `--genome-path-file <GENOME_PATH_FILE>` — File indicating path to genome files to sketch (TSV: genome ID followed by path to FASTA file)
* `-o`, `--output-file <OUTPUT_FILE>` — Output sketch file
* `-w`, `--weighted` — Generated sketch indicating number of times each k-mer occurs
* `-k`, `--kmer-length <KMER_LENGTH>` — Length of k-mers to use

  Default value: `31`
* `-s`, `--scale <SCALE>` — Sketch scaling factor

  Default value: `1000`
* `-t`, `--threads <THREADS>` — Number of threads to use

  Default value: `1`



## `sketchm index`

Create k-mer index from sketches

**Usage:** `sketchm index --sketches <SKETCHES>... --output-file <OUTPUT_FILE>`

###### **Options:**

* `-s`, `--sketches <SKETCHES>` — Genome sketches to index
* `-o`, `--output-file <OUTPUT_FILE>` — Output index file



## `sketchm dist`

Compute distances between genome sketches

**Usage:** `sketchm dist [OPTIONS] --query-sketches <QUERY_SKETCHES>... <--reference-sketches <REFERENCE_SKETCHES>...|--reference-index <REFERENCE_INDEX>>`

###### **Options:**

* `-q`, `--query-sketches <QUERY_SKETCHES>` — Query genome sketches
* `-r`, `--reference-sketches <REFERENCE_SKETCHES>` — Reference genome sketches
* `-i`, `--reference-index <REFERENCE_INDEX>` — Reference k-mer index
* `-o`, `--output-file <OUTPUT_FILE>` — Output file [default: stdout]
* `--min-ani <MIN_ANI>` — Only report ANI values above this threshold [0, 100]

  Default value: `0.01`
* `-t`, `--threads <THREADS>` — Number of threads to use

  Default value: `1`



## `sketchm info`

Display information about sketch file

**Usage:** `sketchm info [OPTIONS] --sketch-file <SKETCH_FILE>`

###### **Options:**

* `-s`, `--sketch-file <SKETCH_FILE>` — Sketch file to query for information
* `-o`, `--output-file <OUTPUT_FILE>` — Output file [default: stdout]



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>

