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
