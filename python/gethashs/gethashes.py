#!/usr/bin/env python3

import os, sys
from os.path import join, getsize
from datetime import datetime

epoch = datetime.utcfromtimestamp(0)

def epoch_seconds(dt):
    return (dt - epoch).total_seconds()

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        argv.append('.')
    argv[0] = ''.join([argv[0].rstrip(os.sep),os.sep])
    print(argv)

    for root, dirs, files in os.walk(argv[0]):
        print(root, "consumes", end=" ")
        print(sum(getsize(join(root, name)) for name in files), end=" ")
        print("bytes in", len(files), "non-directory files")
        if '@eaDir' in dirs:
            dirs.remove('@eaDir')  # don't visit "@eaDir" directories
        relroot = root.rpartition(argv[0])[2]
        #if not relroot: relroot = '.'
        print(root,relroot,dirs)
        #print(files)
        for afile in files:
            rname = os.path.join(root,afile)
            #fname = os.sep.join(filter(None,[relroot,afile]))
            fname = os.path.join(relroot,afile)
            mtime = os.path.getmtime(rname)
            stime = datetime.utcfromtimestamp(mtime).strftime('%Y%m%du%H%M%S')
            rtime = datetime.strptime(''.join([stime,'UTC']),'%Y%m%du%H%M%S%Z')
            print(rname,fname,stime,mtime,epoch_seconds(rtime))

if __name__ == '__main__':
    main()



def isBlank (myString):
    return not (myString and myString.strip())

def isNotBlank (myString):
    return bool(myString and myString.strip())
