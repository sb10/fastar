// Copyright © 2018 Genome Research Limited
// Author: Sendu Bala <sb10@sanger.ac.uk>.
// Using ideas from https://adamdrake.com/faster-command-line-tools-with-go.html
// And based on https://github.com/ekg/freebayes/scripts/coverage_to_regions.py
// by Erik Garrison, Copyright (c) 2010 Erik Garrison, Gabor Marth, released
// under the license below.
//
// This file is part of fastar.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
)

var l = log.New(os.Stderr, "", 0)

func main() {
	// arg handling
	var help = flag.Bool("h", false, "print help text")
	var numRegions, bpPerRegion int64
	flag.Int64Var(&numRegions, "r", 0, "number of regions")
	flag.Int64Var(&bpPerRegion, "b", 0, "number of base pairs per region")
	flag.Parse()

	if *help {
		exitHelp("")
	}

	if len(flag.Args()) != 1 {
		exitHelp("ERROR: you must provide the depth output")
	}

	if numRegions == 0 && bpPerRegion == 0 {
		exitHelp("One of -f or -r must be specified")
	}

	if numRegions != 0 && bpPerRegion != 0 {
		exitHelp("-f and -r are mutually exclusive")
	}

	if numRegions < 0 || bpPerRegion < 0 {
		exitHelp("number must be positive")
	}

	depthPath := flag.Arg(0)

	// calculate bpPerRegion if not supplied
	if bpPerRegion == 0 {
		if depthPath == "-" {
			exitHelp("cannot pipe depth file if -b is not specified")
		}

		// get total coverage (we don't store anything from this pass through
		// the depth file because that would use too much memory)
		totalCoverage := sumColumnThree(depthPath)
		bpPerRegion = int64(math.Ceil(float64(totalCoverage) / float64(numRegions)))
		l.Printf("using %d base pairs per region\n", bpPerRegion)
	}

	// go through the depth file to decide the regions
	calculateRegions(depthPath, bpPerRegion)
}

// exitHelp prints help text and exits 0, unless a message is passed in which
// case it also prints that and exits 1.
func exitHelp(msg string) {
	if msg != "" {
		l.Println(msg)
		fmt.Printf("\n")
	}

	fmt.Printf(`fastar is a fasta region generator.

You use it to generate a list of genomic regions with even sequencing coverage.
This can be used with variant callers to reduce the variance in calling runtime
for each region.

The required command-line argument is the path to the output of 'samtools depth'
on your sequence alignment file. Alternatively this can be '-' to pipe the depth
output in. (Don't use the -a option to 'samtools depth'.)

One of -b or -r is also required. Specifying -r results in it taking almost
twice as long, and doesn't work with piped input, since this results in base
pairs per region being calculated by summing the depths in the depth file,
before going through it again to calculate the regions.

Note also that with -r, the number of regions you get isn't guaranteed. It may
be much higher if you have coverage over many small sequences: the minimum
number of regions will be the number of covered sequences in your reference.

Usage: fastar -[r|b] <int> samtools.depth.output
Options:
  -h          this help text
  -r <int>    desired number of regions
  *or*
  -b <int>    desired total coverage (in base pairs) per region
`)

	if msg != "" {
		os.Exit(1)
	}
	os.Exit(0)
}

func die(err error) {
	l.Printf("ERROR: %s", err.Error())
	os.Exit(1)
}

// mapColsOneAndTwo can be used to parse a reference index (.fai) file, which
// is done in the original freebayes script algorithm, but we don't need it
// here.
// func mapColsOneAndTwo(path string) map[string]int64 {
// 	file, err := os.Open(path)
// 	if err != nil {
// 		die(err)
// 	}
// 	defer func() {
// 		errc := file.Close()
// 		if errc != nil {
// 			l.Printf("WARNING: failed to close %s", path)
// 		}
// 	}()
//
// 	m := make(map[string]int64)
// 	scanner := bufio.NewScanner(file)
// 	for scanner.Scan() {
// 		b := scanner.Bytes()
// 		i := 0
// 		for b[i] != '\t' {
// 			i++
// 		}
// 		key := string(b[0:i])
//
// 		var val int64
// 		for i++; b[i] != '\t'; i++ {
// 			val = val*10 + int64(b[i]) - '0'
// 		}
// 		m[key] = val
// 	}
// 	return m
// }

func sumColumnThree(path string) int64 {
	file, err := os.Open(path)
	if err != nil {
		die(err)
	}
	defer func() {
		errc := file.Close()
		if errc != nil {
			l.Printf("WARNING: failed to close %s", path)
		}
	}()

	var sum int64
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		b := scanner.Bytes()
		i := 0
		for b[i] != '\t' {
			i++
		}
		i++
		for b[i] != '\t' {
			i++
		}
		var depth int
		for i++; i < len(b); i++ {
			depth = depth*10 + int(b[i]) - '0'
		}
		sum += int64(depth)
	}
	return sum
}

func calculateRegions(path string, bpPerRegion int64) {
	var input io.Reader
	if path == "-" {
		input = bufio.NewReader(os.Stdin)
	} else {
		var err error
		input, err = os.Open(path)
		if err != nil {
			die(err)
		}
		defer func() {
			errc := input.(*os.File).Close()
			if errc != nil {
				l.Printf("WARNING: failed to close %s", path)
			}
		}()
	}

	var prevSeq string
	var startPos, lastPos, bpInRegion int64

	scanner := bufio.NewScanner(input)
	for scanner.Scan() {
		// parse the line
		b := scanner.Bytes()
		i := 0
		for b[i] != '\t' {
			i++
		}
		seq := string(b[0:i])

		var pos int64
		for i++; b[i] != '\t'; i++ {
			pos = pos*10 + int64(b[i]) - '0'
		}

		var depth int
		for i++; i < len(b); i++ {
			depth = depth*10 + int(b[i]) - '0'
		}

		// print regions that end on sequence changes or when over bpPerRegion.
		// region coordinates are 1-based, inclusive of start and end.
		if prevSeq != seq {
			if prevSeq == "" {
				prevSeq = seq
				startPos = pos
			} else {
				fmt.Printf("%s:%d-%d\n", prevSeq, startPos, lastPos)
				startPos = pos
				prevSeq = seq
				bpInRegion = 0
			}
		}

		bpInRegion += int64(depth)
		if bpInRegion > bpPerRegion {
			fmt.Printf("%s:%d-%d\n", seq, startPos, pos)
			startPos = pos + 1
			bpInRegion = 0
		}

		lastPos = pos
	}

	if lastPos > startPos {
		fmt.Printf("%s:%d-%d\n", prevSeq, startPos, lastPos)
	}
}
