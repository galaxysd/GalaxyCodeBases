#!/usr/bin/env python3
import sys

def main():
    if len(sys.argv) == 1 :
        print('[x]No arguments given!',file=sys.stderr,flush=True)

        print("Name of the script:",sys.argv[0])
        print("Number of arguments:", len(sys.argv))
        print("Arguments:" , str(sys.argv))
    if len(sys.argv) < 3 :
        print('Usage:',sys.argv[0],'<in.bam> <outprefix> [verbose=0]',file=sys.stderr,flush=True)
        exit(0)
    try:
        verbose = int(sys.argv[3])
    except ValueError:
        verbose = 0
    except IndexError:
        verbose = 0

    inBAMname = sys.argv[1]
    outPrefix = sys.argv[2]
    print('From:[{}], To:[{}.(Watson|Crick|Dropped).bam].\nVerbose: [{}].'.format(inBAMname,outPrefix,verbose),file=sys.stderr,flush=True)
    splitBSbam(inBAMname,outPrefix,3,verbose)

class DebugWrite:
    def __init__(self, name=None, verbose=1):
        self._name = name
        self.verbose = verbose

    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, name):
        self._name = name

    def write(self, read):
        if self.verbose > 1:
            print ("[{}]\t{}\t{}".format(self.name, read.query_name, read.flag))
    def close(self):
        return 0

def splitBSbam(inBAMname,outPrefix,n_threads=3,verbose=0):
    import pysam
    samfile = pysam.AlignmentFile(inBAMname, "rb")
    if verbose == 0:
        fileWatson = pysam.AlignmentFile('.'.join((outPrefix,'Watson','bam')), "wb", template=samfile, threads=n_threads)
        fileCrick = pysam.AlignmentFile('.'.join((outPrefix,'Crick','bam')), "wb", template=samfile, threads=n_threads)
        fileDropped = pysam.AlignmentFile('.'.join((outPrefix,'Dropped','bam')), "wb", template=samfile, threads=n_threads)
    else:
        fileWatson = DebugWrite('oWatson',verbose)
        fileCrick = DebugWrite('oCrick',verbose)
        fileDropped = DebugWrite('oDropped',verbose)
    from collections import defaultdict
    stats = defaultdict(int)
    stats['Watson']=stats['Crick']=0    # make them be first in `str(stats)`.
    for read in samfile:
        if read.flag & 0xF00 :
            fileDropped.write(read)
            stats['Dropped_Supplementary'] += 1
            continue
        try:
            tagYD = read.get_tag('YD')
            stats['etc_YD'] += 1
            if tagYD == 'f':
                fileWatson.write(read)
                stats['Watson'] += 1
            elif tagYD == 'r':
                fileCrick.write(read)
                stats['Crick'] += 1
            else:
                print(tagYD,read.query_name,read.flag)
                exit(3)
            #continue
        except KeyError:
            #print('[!]',sys.exc_info()[1])    # (<class 'KeyError'>, KeyError("tag 'YD' not present"), <traceback object at 0x10f664ec8>)
            #pass
            if read.is_unmapped:
                stats['etc_nonYD'] += 1
                if read.mate_is_unmapped:
                    fileDropped.write(read)
                    stats['Dropped_PEUnmapped'] += 1
                elif (read.is_read2 and (not read.mate_is_reverse)) or (read.is_read1 and read.mate_is_reverse) :
                    fileWatson.write(read)
                    stats['Watson'] += 1
                elif (read.is_read2 and read.mate_is_reverse) or (read.is_read1 and (not read.mate_is_reverse)) :
                    fileCrick.write(read)
                    stats['Crick'] += 1
                else:
                    print (str(read))
                    exit(4)
            else:   # not bwa-meth
                stats['etc_nonBWA-meth'] += 1
                if (read.is_read1 and (not read.is_reverse)) or (read.is_read2 and read.is_reverse) :
                    fileWatson.write(read)
                    stats['Watson'] += 1
                elif (read.is_read1 and read.is_reverse) or (read.is_read2 and (not read.is_reverse)) :
                    fileCrick.write(read)
                    stats['Crick'] += 1
                else:
                    print (str(read))
                    exit(2)
    else:
        print(str(stats),file=sys.stderr,flush=True)
        #import pprint
        #pprint.pprint(stats, indent=2, stream=sys.stderr)
    fileDropped.close()
    fileWatson.close()
    fileCrick.close()
    samfile.close()

if __name__ == "__main__":
    main()  # time ./samflagfilter.py t.bam tt2
