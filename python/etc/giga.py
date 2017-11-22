#!/usr/bin/env python3

import os, sys, bz2, csv, re, json, sqlite3
import urllib.request

import pprint

myUA='Mozilla/5.0 (Macintosh; Intel Mac OS X 10.12; rv:58.0) Gecko/20100101 Firefox/58.0'
conn = sqlite3.connect('giga.authors.sqlite')
conn.execute('CREATE TABLE PubDat (DOI TEXT, Title TEXT, Type TEXT, Authors TEXT, RefList TEXT)')


with bz2.open('giga.tsv.bz2','rt',encoding='utf-8') as tsvin:
    tsvin = csv.reader(tsvin, delimiter='\t')
    for row in tsvin:
        if row[0] == 'Title' or len(row)==0: continue
        print((len(row),row[0]))
        theurl = 'https://academic.oup.com/gigascience/article-lookup/doi/' + row[16];
        req = urllib.request.Request(theurl)
        req.add_header('Referer', 'https://academic.oup.com/gigascience/')
        req.add_header('User-Agent',myUA)
        with urllib.request.urlopen(req) as r:
            htm = r.read().decode('utf-8')
        it = iter( htm.split('\n') )
        data={'strAuthors':'NA','reflist':'NA','tocSections':'NA'}
        for line, secline in zip(it, it):
            if re.search('<script type="application\/ld\+json">',line):
                datAuthors = json.loads(secline)
                data['strAuthors'] = json.dumps(datAuthors['author'])
                #pprint.pprint(datAuthors['author'])
            if re.search('<div class="ref-list">',line):
                data['reflist'] = line.strip('\t\r \n')
            if re.search('Issue Section',line):
                m = re.search('>([^<>]+)<\/a>',secline)
                if m:
                    data['tocSections'] = m.group(1)
        #pprint.pprint(data)
        conn.execute('INSERT INTO PubDat ( DOI,Title,Type,Authors,RefList ) VALUES ( ?,?,?,?,? )',(row[16],row[0], data['tocSections'],data['strAuthors'],data['reflist'] ) )
        conn.commit()

conn.close()
