#!/usr/bin/env python3

import os, sys, hashlib
import sqlite3, getopt, mmap
#from os.path import join, getsize
from datetime import datetime

epoch = datetime.utcfromtimestamp(0)
def epoch_seconds(dt):
    return (dt - epoch).total_seconds()

BUF_SIZE = 1048576  # lets read stuff in 1Mb chunks!
def sha1file(fname=None, blocksize=BUF_SIZE):
    if fname is None:
        return None
    sha1 = hashlib.sha1()
    with open(fname, 'rb') as f:
        with mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ) as mm:
            for block in iter(lambda: mm.read(blocksize), b""):
                sha1.update(block)
    return sha1.hexdigest()

class Config: # https://stackoverflow.com/a/47016739/159695
    def __init__(self, **kwds):
        self.verbose=0 # -1=quiet  0=norm  1=noisy
        self.mode=0 # 
        self.output = ''.join(['.',os.sep])
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
    phelp('Usage: gethashes [opts] [-p dir] [-T|-C] [-f file] [files...]')
    phelp('  -T       Test mode (default)')
    phelp('  -C       Create mode')
    #phelp('  -t <t>   set type to <t> (%s, or auto(default))'%', '.join(sorted(hashlib.algorithms_available)))
    phelp('  -f <f>   use <f> as list file (.hash)')
    phelp('  -p <d>   change to directory <d> before doing anything')
    phelp('Options in Create mode:')
    phelp('  -u [s][k,m]   load .sha1 files in subdirectories and skip older recorded files larger than [s] [*1024, *1048576] (default=1m)')
    phelp('  -s            Always skip recorded files even if loaded .sha1 file is older')
    phelp('  -1            Also create .sha1 file')
    phelp('Options in Test mode:')
    phelp('  -b <f>        Output list of bad files to file <f>')
    phelp('Other Options:')
    phelp('  -v/-q    verbose/quiet, change verbosity [-1,2]')
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
        opts, args = getopt.gnu_getopt(argv, "CTf:p:u:s1b:vqh?", ['help','version'])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        printusage(1)

    try:
        prevopt=''
        for o, a in opts:
            if o=='-p':
                config.startpoint = ''.join([a.rstrip(os.sep),os.sep])
                os.chdir(a)
            elif o in ("-o", "--output"):
                config.output = a
            elif o=='-v':
                if config.verbose >=0: config.verbose +=1
            elif o=='-q':
                if config.verbose >=0:
                    config.verbose =-1
                else:
                    config.verbose -=1
            elif o=='-h' or o=='-?' or o=='--help':
                printusage()
            elif o=='-V' or o=='--version':
                print('gethashes %s'%version)
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
    print(config)

    for root, dirs, files in os.walk(config.output): # os.walk(top, topdown=True, onerror=None, followlinks=False)
        #print(root, "consumes", end=" ")
        #print(sum(getsize(join(root, name)) for name in files), end=" ")
        #print("bytes in", len(files), "non-directory files")
        if '@eaDir' in dirs:
            dirs.remove('@eaDir')  # don't visit "@eaDir" directories
        print(dirs,files)
        dirs.sort(reverse=True)
        files.sort(reverse=True)
        print(dirs,files)
        relroot = root.rpartition(config.output)[2]
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
