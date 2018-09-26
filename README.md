# fastar

fastar is a fasta region generator.

You use it to generate a list of genomic regions with even sequencing coverage.
This can be used with variant callers to reduce the variance in calling runtime
for each region.

# Installation

The easiest way to install fastar is to download the latest binary from
[here](https://github.com/sb10/fastar/releases) and then `chmod +x` the file.

To build from source you need a working go installation. Then just:

```
go install github.com/sb10/fastar
```

# Usage

The command has help accessible by running with the `-h` option.

```
fastar -[r|b] <int> samtools.depth.output
```

It is recommended to use -b to set the desired total coverage (in base pairs)
per region, eg. 2000000. Using -r instead results in it taking almost twice as
long, and doesn't work with piped input, since this results in base pairs per
region being calculated by summing the depths in the depth file, before going
through it again to calculate the regions.

Note also that with -r, the number of regions you get isn't guaranteed. It may
be much higher if you have coverage over many small sequences: the minimum
number of regions will be the number of covered sequences in your reference.

# Provenance

fastar is a pure-go version of the 'coverage_to_regions.py' script supplied with
[freebayes](https://github.com/ekg/freebayes).

fastar was written because the python version is slow and uses lots of memory:
on this author's ~3billion line depth file, it consumed over 250GB of memory
before crashing after 23hrs.

fastar, on the same file, takes ~12mins (~8mins if using -b) and consumes a
constant small amount of memory. On a ~50million line truncated version of the
file, coverage_to_regions.py took 308s (and used ~1GB), while fastar took 14s
(and used 7MB), 9s with -b, making it about 22-34 times faster.

fastar differs from coverage_to_regions.py in that its regions start and end
on the first and last base per reference sequence that had coverage, instead of
starting at 0 (an invalid coordinate) and ending at the end of the sequence.
Because of this you don't need to supply the reference index file. Regions
output also do not overlap.

There is a tool in [goleft](https://github.com/brentp/goleft) called indexsplit
that operates in seconds and doesn't require the depth file, but it isn't as
accurate. However, typically accuracy isn't that important for this use case,
and indexsplit may be good enough.

That said, what determines calling time is related to base pairs per region, so
it is more useful to generate regions based on a fixed value for that, as
opposed to a fixed number of regions, which is the only option available with
goleft (and the original freebayes script).

fastar also only outputs regions that have coverage, so you don't get regions
for whole references sequences where no alignments occurred. This makes the
output more suitable for passing to callers, since you now won't end up wasting
time running calls on zero data, and you don't have to deal with empty VCF
files.
