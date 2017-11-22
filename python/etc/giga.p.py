#!/usr/bin/env python3

import os, sys, bz2, csv, re, json, sqlite3

conn = sqlite3.connect('giga.authors.sqlite')

for row in conn.execute("SELECT * FROM PubDat"):
    #print(row)
    auther=json.loads(row[3])
    print('DOI:%s\nTitle:%s\nType:%s\nRef:{%s}\nAuthers:'%(row[0],row[1],row[2],row[4]))
    print(auther)

