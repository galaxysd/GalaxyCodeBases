#!/usr/bin/env python3
import pprint

def main():
    import sys
    print("This is the name of the script: ",sys.argv[0])
    print("Number of arguments: ", len(sys.argv))
    print("The arguments are: " , str(sys.argv))
    if len(sys.argv) == 1 :
        print('[x]No arguments given!',file=sys.stderr,flush=True)
        exit(0)
    print('Reading',sys.argv[1])
    inBAMname = sys.argv[1]
    import pysam
    samfile = pysam.AlignmentFile(inBAMname, "rb")
    pairedreads = pysam.AlignmentFile("allpaired.bam", "wb", template=samfile)
    for read in samfile:
         if read.is_paired:
                 #pairedreads.write(read)
                 print (str(read))
                 pprint.pprint(read)

    pairedreads.close()
    samfile.close()

if __name__ == "__main__":
    main()
