#!/usr/bin/env python3

import os, sys, hashlib
import sqlite3, getopt
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

class Config:
    def __init__(self, **kwds):
        self.verbose=0 # -1=quiet  0=norm  1=noisy
        self.__dict__.update(kwds) # Must be last to accept assigned member variable.
    def __repr__(self):
        args = ['%s=%s' % (k, repr(v)) for (k,v) in vars(self).items()]
        return '%s(%s)' % ( self.__class__.__qualname__, ', '.join(args) )

def pverbose(s,nl='\n'):
    if config.verbose>0:
        sys.stdout.write(s+nl)
def pinfo(s,nl='\n'):
    if config.verbose>=0 or config.verbose==-3:
        sys.stdout.write(s+nl)
def perror(s,nl='\n'):
    if config.verbose>=-1:
        sys.stdout.flush() # avoid inconsistent screen state if stdout has unflushed data
        sys.stderr.write(s+nl)
def printusage(err=0):
    phelp = err and perror or pinfo # False->pinfo, True->perror
    phelp('Usage: gethashes [opts] [-p dir] [-T|-C] [-t type] [-f file] [files...]')
    phelp('  -p <d>   change to directory <d> before doing anything')
    phelp(' --help/-h show help')
    phelp(' --version show gethashes and module versions')
    sys.exit(err)

config=Config()

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        argv.append('.')
        argv[0] = ''.join([argv[0].rstrip(os.sep),os.sep])
    print(argv)

    try:
        opts, args = getopt.gnu_getopt(argv, "ho:v?", ['help','version', "output="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        printusage(1)

    try:
        prevopt=''
        for o, a in opts:
            if o=='-p':
                chdir(a)
            elif o in ("-o", "--output"):
                output = a
            elif o=='-v':
                config.setx('verbose', 'v')
            elif o=='-h' or o=='-?' or o=='--help':
                printusage()
            elif o=='--version':
                print('cfv %s'%version)
                try:
                    if not nommap: print('+mmap')
                except NameError: pass
                try: print('fchksum %s'%fchksum.version())
                except NameError: pass
                print('python %08x-%s'%(sys.hexversion,sys.platform))
                sys.exit(0)
            else:
                assert False, "unhandled option"
            prevopt=o
    except RuntimeError as e:
        perror('cfv: %s'%e)
        sys.exit(1)

    for root, dirs, files in os.walk(argv[0]): # os.walk(top, topdown=True, onerror=None, followlinks=False)
        #print(root, "consumes", end=" ")
        #print(sum(getsize(join(root, name)) for name in files), end=" ")
        #print("bytes in", len(files), "non-directory files")
        if '@eaDir' in dirs:
            dirs.remove('@eaDir')  # don't visit "@eaDir" directories
        print(dirs,files)
        dirs.sort(reverse=True)
        files.sort(reverse=True)
        print(dirs,files)
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

version='0.1'
if __name__ == '__main__':
    main()

def isBlank (myString):
    return not (myString and myString.strip())

def isNotBlank (myString):
    return bool(myString and myString.strip())
