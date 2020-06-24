#!/usr/bin/env python

class PhredHelper:
    def __init__(self, offset=33):
        if offset != 33 and offset != 64:
            sys.err.write("PhredHelper: invalid offset specified. Using default (Sanger 'Q+33')")
            self.offset = 33
        else:
            self.offset = offset

    def char_to_int(self, char):
        raw_val = ord(char)
        phred_val = raw_val - self.offset
        return phred_val

    def int_to_char(self, num):
        return chr(num + self.offset)
