#!/usr/bin/env python3

import os, sys, hashlib
import sqlite3, getopt, mmap
import re
import gzip
#from os.path import join, getsize
from datetime import datetime
import pprint

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
        self.mode=0 # 0:Test, 1:Create
        self.skipNewer=0
        self.sha1dump=0
        self.ssize=1048576
        self.startpoint = ''.join(['.',os.sep])
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
    #phelp('  -T       Test mode (default)')
    #phelp('  -C       Create mode')
    #phelp('  -t <t>   set type to <t> (%s, or auto(default))'%', '.join(sorted(hashlib.algorithms_available)))
    phelp('  -p <d>   change to directory <d> before doing anything')
    phelp('  -f <f>   use <f> as list file (<d>.hash.gz)')
    #phelp('Options in Create mode:')
    phelp('  -s [s][k,m]   load .sha1 files in subdirectories and skip older recorded files larger than [s] [*1024, *1048576] (default=1m)')
    phelp('  -a            Always skip recorded files even if loaded .sha1 file is older')
    #phelp('  -1            Also create <f>.sha1 file')
    #phelp('Options in Test mode:')
    #phelp('  -b <l>        Output list of bad files to file <l>')
    #phelp('Other Options:')
    phelp('  -v/-q    verbose/quiet, change verbosity [-1,2]')
    phelp(' --help/-h show help')
    phelp(' --version show gethashes and module versions')
    sys.exit(err)

config=Config()

# https://stackoverflow.com/questions/635483/what-is-the-best-way-to-implement-nested-dictionaries/19829714#19829714
class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)() # retain local pointer to value
        return value                     # faster to return than dict lookup

OldHashes = Vividict()
HitHashes = 0
# rem from `cfv` L1116:`_foosum_rem`. `re.match()` checks for a match only at the beginning of the string, thus not r'^'.
sha1rem=re.compile(r'([0-9a-fA-F]{40}) ([ *])([^\r\n]+)[\r\n]*$')

def loadsha1(root,afile):
    global OldHashes
    rname = os.path.join(root,afile)
    mtime = os.path.getmtime(rname)
    itextmode = 0
    for line in open(rname):
        x = sha1rem.match(line)
        if not x: return -1
        if x.group(2)==' ':
            if not itextmode:
                pinfo('[!]Textmode in "%s".'%(rname))
            itextmode += 1
            continue
        #pprint.pprint([x.group(3),x.group(1),x.group(2)])
        iname = os.path.join(root,x.group(3))
        istat = os.stat(iname)
        if istat.st_size < config.ssize:
            continue
        if not config.skipNewer:
            itime = os.path.getmtime(iname)
            #isize = os.path.getsize(iname)
            pprint.pprint(['t:',iname,mtime,itime])
            #pprint.pprint(['->',iname,itime,istat.st_size,afile,istat.st_ino,istat.st_dev])
            if mtime < itime:
                continue
        OldHashes[istat.st_dev][istat.st_ino] = x.group(1)
        #pprint.pprint(('O:',OldHashes))
    if itextmode>1 :
        pinfo('[!]Textmode %d times in "%s" !'%(itextmode,rname))
    return

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        argv.append('.')
        argv[0] = ''.join([argv[0].rstrip(os.sep),os.sep])
    pprint.pprint(argv) # <-- DEBUG

    try:
        opts, args = getopt.gnu_getopt(argv, "CTf:p:s:a1b:vqh?", ['help','version'])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err)  # will print something like "option -a not recognized"
        printusage(1)

    try:
        prevopt=''
        for o, a in opts:
            if o=='-p':
                config.startpoint = ''.join([a.rstrip(os.sep),os.sep])
                os.chdir(a) # also checks PermissionError and FileNotFoundError for me
            elif o=='-f':
                config.hashfile = a
            elif o=='-s':
                config.ssize = human2bytes(a)
            elif o=='-a':
                config.skipNewer = 1
            elif o=='-1':
                config.sha1dump = 1
            elif o=='-v':
                if config.verbose >=0: config.verbose +=1
            elif o=='-q':
                if config.verbose >=0:
                    config.verbose =-1
                else:
                    config.verbose -=1
            elif o in ("-h", "--help", '-?'):
                printusage()
            elif o=='-V' or o=='--version':
                print('gethashes %s'%version)
                print('python %08x-%s'%(sys.hexversion,sys.platform))
                sys.exit(0)
            else:
                assert False, "unhandled option"
            prevopt=o
    except RuntimeError as e:
        perror('cfv: %s'%e)
        sys.exit(1)
    if not hasattr(config, 'hashfile'):
        dirName = os.path.basename(os.path.abspath(config.startpoint))
        #dirName = os.path.basename(os.getcwd())
        config.hashfile = ''.join([config.startpoint, dirName, '.hash.gz'])
    pprint.pprint(config) # <-- DEBUG
    doCreation()

def doCreation():
    f_out = gzip.open(config.hashfile, 'wt', encoding='utf-8')

    for root, dirs, files in os.walk(config.startpoint): # os.walk(top, topdown=True, onerror=None, followlinks=False)
        #print(root, "consumes", end=" ")
        #print(sum(getsize(join(root, name)) for name in files), end=" ")
        #print("bytes in", len(files), "non-directory files")
        if '@eaDir' in dirs:
            dirs.remove('@eaDir')  # don't visit "@eaDir" directories
        #print(dirs,files)
        dirs.sort(reverse=True)
        files.sort(reverse=True)
        #print(dirs,files)
        relroot = root.rpartition(config.startpoint)[2]
        #if not relroot: relroot = '.'
        #print(root,relroot,dirs)
        #print(files)
        global OldHashes
        global HitHashes
        for afile in files:
            if afile.endswith(".sha1"):
                loadsha1(root,afile)
        for afile in files:
            rname = os.path.join(root,afile)
            if os.path.samefile(rname,config.hashfile):
                continue
            #fname = os.sep.join(filter(None,[relroot,afile]))
            fname = os.path.join(relroot,afile)
            istat = os.stat(rname)
            if (istat.st_dev in OldHashes) and (istat.st_ino in OldHashes[istat.st_dev]):
                ihash = OldHashes[istat.st_dev][istat.st_ino]
                HitHashes += 1
            else:
                ihash = sha1file(rname)
            pprint.pprint(('O:',fname,rname,ihash,HitHashes))
            f_out.write('%s *%s\n'%(ihash,fname))
            #mtime = os.path.getmtime(rname)
            #stime = datetime.utcfromtimestamp(mtime).strftime('%Y%m%du%H%M%S')
            #rtime = datetime.strptime(''.join([stime,'UTC']),'%Y%m%du%H%M%S%Z')
            #print(rname,fname,stime,mtime,epoch_seconds(rtime))
            #fsize = os.path.getsize(rname)
            #hsize = bytes2human(fsize)
            #print(fname,hsize,stime,sha1file(rname),sep='\t')
    f_out.close()
    if HitHashes:
        pinfo('[!]Skipped hashing of %d recorded file(s).'%HitHashes)
    pinfo('[!]Done. Test with `cfv -T -f %s`.'%config.hashfile)

# https://github.com/giampaolo/pyftpdlib/blob/0430c92e9d852a6d175b489c0ebf17fbc0190914/scripts/ftpbench#L139
def bytes2human(n, format="%(value).1f%(symbol)s", intfmt="%(value).0f %(symbol)s"):
    """
    >>> bytes2human(10000)
    '9K'
    >>> bytes2human(100001221)
    '95M'
    """
    symbols = ('B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols[1:]):
        prefix[s] = 1 << (i + 1) * 10
    for symbol in reversed(symbols[1:]):
        if n >= prefix[symbol]:
            value = float(n) / prefix[symbol]
            return format % locals()
    #import re
    #intfmt=re.sub(r'\(value\)\.(\d+)',r'(value).0',format)
    #print(intfmt)
    return intfmt % dict(symbol=symbols[0], value=n)

# http://goo.gl/zeJZl
def human2bytes(s):
    """
    >>> human2bytes('1M')
    1048576
    >>> human2bytes('1G')
    1073741824
    """
    symbols = ('B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    letter = s[-1:].strip().upper()
    num = s[:-1]
    assert num.isdigit() and letter in symbols, s
    num = float(num)
    prefix = {symbols[0]: 1}
    for i, s in enumerate(symbols[1:]):
        prefix[s] = 1 << (i + 1) * 10
    return int(num * prefix[letter])

version='0.1'
if __name__ == '__main__':
    main()

def isBlank (myString):
    return not (myString and myString.strip())

def isNotBlank (myString):
    return bool(myString and myString.strip())
