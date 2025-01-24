# Copyright 2016 Kevin Murray <spam@kdmurray.id.au>
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from docopt import docopt
from pymer import CountMinKmerCounter
import screed
from sys import stderr, stdout
from time import time

def benchmark_counting():
    cli = '''

    USAGE:
        benchmark.py [options] OUTFILE READFILES ...

    OPTIONS:
        -k KSIZE    Kmer length [default: 20]
        -N NTAB     Number of tables [default: 4]
        -x TSIZE    Table size [default: 1e9]
    '''

    opts = docopt(cli)
    k = int(opts['-k'])
    nt = int(opts['-N'])
    ts = int(float(opts['-x']))
    outfile = opts['OUTFILE']
    readfiles = opts['READFILES']

    counter = CountMinKmerCounter(k, sketchshape=(nt, ts))

    start = time()
    nread = 0
    for readfile in readfiles:
        print("Consuming:",  readfile)
        with screed.open(readfile) as reads:
            for i, read in enumerate(reads, start=1):
                counter.consume(str(read.sequence))
            nread += i
    endhash = time()
    took = endhash - start
    print('Counted', nread, 'reads in {:0.1f} seconds'.format(took),
          '({:0.2f} K reads/sec)'.format((nread / 1000)/took))

    print("Writing to", outfile)
    counter.write(outfile)
    print('Done, writing took {:0.1f} sec'.format(time() - endhash))


if __name__ == '__main__':
    benchmark_counting()
