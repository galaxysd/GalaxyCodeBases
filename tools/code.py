#!/usr/bin/python
import sys
import re
# Format
# @I321_1_FC30BBUAAXX:2:1:4:1436/1
ver=0
code = {}
for i in xrange(0,10):
    code[str(i)] = i
for i in xrange(65,91):
    code[chr(i)] = i-55
for i in xrange(97,123):
    code[chr(i)] = i-61
rev_code = [0,]*62
for (key, value) in code.iteritems():
    rev_code[value]=key
for i in xrange(62):
    print i,rev_code[i]

def int2bin(integer):
    bin = ""
    for i in xrange(63,-1,-1):
        bin += str((integer>>i)&1)
    return bin
#regexpNew = re.compile("@FC(\w{1,5})AAXX:(\d):(\d+):(-?\d+):(-?\d+)#[01]/?([12]?)")
regexpOld = re.compile("@(\w+)_FC(\w{1,5})AAXX:(\d):(\d+):(-?\d+):(-?\d+)/?([12]?)")
#regexpNew = re.compile("@FC(\w{1,5})AAXX:(\d):(\d+):(-?\d+):(-?\d+)#[01]/?([12]?)")
regexpNew = re.compile("@FC(\w{1,5})AAXX:(\d):(\d+):(-?\d+):(-?\d+)#\w+/?([12]?)")
def scode(line):
    try:
        data = ''
        if re.match(regexpNew, line).groups():
           data=re.match(regexpNew, line).groups()
        elif re.match(regexpOld, line).groups():
           data=re.match(regexpOld, line).groups()
    except AttributeError:
        print line.rstrip()
        sys.exit(255)
    pe = '1'
    if re.match(regexpNew, line).groups():
        if type(data) == None:
            print >>sys.stderr, "The FASTQ header cannot be recognized\n",line
            sys.exit(255)
        elif len(data) == 5:
            ( fc, lane, tile, x, y) = data
        elif len(data) == 6:
            ( fc, lane, tile, x, y, pe) = data
        else:
            print >>sys.stderr, "The FASTQ header cannot be recognized\n",line
            sys.exit(255)
    elif re.match(regexpOld, line).groups():
        if type(data) == None:
            print >>sys.stderr, "The FASTQ header cannot be recognized\n",line
            sys.exit(255)
        elif len(data) == 6:
            (machine, fc, lane, tile, x, y) = data
        elif len(data) == 7:
            (machine, fc, lane, tile, x, y, pe) = data
        else:
            print >>sys.stderr, "The FASTQ header cannot be recognized\n",line
            sys.exit(255)
    fc_num = 0
    for char in fc:
        fc_num = (fc_num<<6)|code[char]
    #print fc_num,int2bin(fc_num)
    (lane, tile, x, y) = map(int, (lane, tile, x, y))
    if x < 0:
        x += 4096
    if y < 0:
        y += 4096
    #print "FC"+fc+"AAXX",lane,tile,x,y,pe
    lane-=1; tile-=1;#pe-=1
    assert lane<8 and tile<512 and x < 4096 and y < 4096 #and pe<2
    info = (lane<<34)|(tile<<25)|(x<<13)|(y<<1)|1
    #print info,int2bin(info)
    part1 = (fc_num<<34)|(info>>3&(~(~0<<34)))
    #print part1,int2bin(part1)
    part2 = (ver<<3)|(info&7)
    #print part2,int2bin(part2)
    code_string = "@"
    while part2 != 0:
        code_string += rev_code[part2 % 62]
        part2 /= 62
        #print part2,code_string
    code_string += "_"
    while part1 != 0:
        code_string += rev_code[part1 % 62]
        part1 /= 62
        #print part1,code_string
    code_string += ("/"+pe)
    return code_string

def decode(code_string):
    part1 = part2 = 0
    is_part1 = True
    for char in code_string[::-1]:
        if char == "_":
            is_part1 = False
        elif is_part1:
            part1 = part1*62+code[char]
        else:
            part2 = part2*62+code[char]
    #print part1,part2
    fc_num = (part1>>34)&(~(~0<<30))
    #print fc_num
    lane = ((part1>>31)&7)+1
    tile = ((part1>>22)&511)+1
    x = (part1>>10)&4095
    y = (part1&1023)<<2|((part2>>1)&3)
    pe = (part2&1)+1
    fc = ""
    while fc_num != 0:
        fc = rev_code[fc_num&63] + fc
        fc_num >>= 6
    #print "FC"+fc+"AAXX",lane,tile,x,y,pe

try:
    fastq = open(sys.argv[1])
    output = open(sys.argv[2],"w")
except IndexError:
    print >>sys.stderr, "python",sys.argv[0],"FASTQ OUTPUT"
    sys.exit(1)
except IOError:
    print >>sys.stderr, "No such file or directory:",sys.argv[1]
    sys.exit(1)

line_count = 0
for line in fastq:
    line_count += 1
    if line_count % 4 == 1:
        temp = scode(line)
        print >>output, temp
        #decode(temp)
    elif line_count % 4 == 3:
        print >>output, "+"
    else:
        print >>output, line.rstrip()
