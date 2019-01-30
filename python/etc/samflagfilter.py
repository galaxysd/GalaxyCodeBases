#!/usr/bin/env python3

def main(n_threads):
    import sys
    if len(sys.argv) == 1 :
        print('[x]No arguments given!',file=sys.stderr,flush=True)

        print("Name of the script:",sys.argv[0])
        print("Number of arguments:", len(sys.argv))
        print("Arguments:" , str(sys.argv))
    if len(sys.argv) != 3 :
        print('Usage:',sys.argv[0],'<in.bam> <outprefix>',file=sys.stderr,flush=True)
        exit(0)
    inBAMname = sys.argv[1]
    outPrefix = sys.argv[2]
    print('Read From:',inBAMname)
    print('Write To:','.'.join((outPrefix,'(Watson|Crick|Supplementary)','bam')))

    import pysam
    samfile = pysam.AlignmentFile(inBAMname, "rb")
    fileWatson = pysam.AlignmentFile('.'.join((outPrefix,'Watson','bam')), "wb", template=samfile, threads=n_threads)
    fileCrick = pysam.AlignmentFile('.'.join((outPrefix,'Crick','bam')), "wb", template=samfile, threads=n_threads)
    fileSupplementary = pysam.AlignmentFile('.'.join((outPrefix,'Supplementary','bam')), "wb", template=samfile, threads=n_threads)
    for read in samfile:
        if read.flag & 0xF00 :
            fileSupplementary.write(read)
        elif ((read.flag & 0x0040) and (not (read.flag & 0x0010))) or ((read.flag & 0x0080) and (read.flag & 0x0010)) :
            fileWatson.write(read)
        elif ((read.flag & 0x0040) and (read.flag & 0x0010)) or ((read.flag & 0x0080) and (not (read.flag & 0x0010))) :
            fileCrick.write(read)
        else:
            print (str(read))
            exit(2)
    fileSupplementary.close()
    fileWatson.close()
    fileCrick.close()
    samfile.close()

if __name__ == "__main__":
    main(3)
