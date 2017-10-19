#!/usr/bin/env python3

import os, sys, hashlib
#from os.path import join, getsize
from datetime import datetime

epoch = datetime.utcfromtimestamp(0)
def epoch_seconds(dt):
    return (dt - epoch).total_seconds()

BUF_SIZE = 1048576  # lets read stuff in 1Mb chunks!
def sha1file(fname=None):
    if fname is None:
        return None
    sha1 = hashlib.sha1()
    with open(fname, 'rb') as f:
        while True:
            data = f.read(BUF_SIZE)
            if not data:
                break
            sha1.update(data)
    return sha1.hexdigest()

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        argv.append('.')
    argv[0] = ''.join([argv[0].rstrip(os.sep),os.sep])
    #print(argv)

    for root, dirs, files in os.walk(argv[0]):
        #print(root, "consumes", end=" ")
        #print(sum(getsize(join(root, name)) for name in files), end=" ")
        #print("bytes in", len(files), "non-directory files")
        if '@eaDir' in dirs:
            dirs.remove('@eaDir')  # don't visit "@eaDir" directories
        relroot = root.rpartition(argv[0])[2]
        #if not relroot: relroot = '.'
        #print(root,relroot,dirs)
        #print(files)
        for afile in files:
            rname = os.path.join(root,afile)
            #fname = os.sep.join(filter(None,[relroot,afile]))
            fname = os.path.join(relroot,afile)
            mtime = os.path.getmtime(rname)
            stime = datetime.utcfromtimestamp(mtime).strftime('%Y%m%du%H%M%S')
            #rtime = datetime.strptime(''.join([stime,'UTC']),'%Y%m%du%H%M%S%Z')
            #print(rname,fname,stime,mtime,epoch_seconds(rtime))
            print(fname,stime,sha1file(rname),sep='\t')

if __name__ == '__main__':
    main()



def isBlank (myString):
    return not (myString and myString.strip())

def isNotBlank (myString):
    return bool(myString and myString.strip())
