#!/usr/bin/env python3
import sys
import css_parser
import re
from pprint import pprint

CSSstr='''
@font-face {
	font-family: "Arial";
	src: url("https://www.com/test1.eot");
	src: url("https://www.com/test1.eot?#iefix") format("embedded-opentype"), url("https://www.com/test1.woff2") format("woff2"), url("https://www.com/test1.woff") format("woff"), url("https://www.com/test1.svg") format("svg");
}
@font-face {
	font-family: "Arial Italic";
	src: url("www.com/test2.eot");
	src: url("www.com/test2.eot?#iefix") format("embedded-opentype"), url("www.com/test2.woff2") format("woff2"), url("www.com/test2.woff") format("woff"), url("www.com/test2.svg") format("svg");
}
'''

parser = css_parser.CSSParser()
sheet = parser.parseFile('WF-000099-001641.css')
#sheet = parser.parseString(CSSstr)

fonts = {}

for rule in sheet:
    if rule.type == rule.FONT_FACE_RULE:
        # find property
        for property in rule.style:
            if property.name == 'font-family':
                font = property.value.strip('"')
                fonts[font] = {}
            if property.name == 'src':
                for i in range(0, property.propertyValue.length, 2):
                    url = property.propertyValue[i].absoluteUri
                    format = re.findall('^format\("(.*?)"\)$', property.propertyValue[i + 1].value)[0]
                    #fonts[font][pair[1]] = pair[0]
                    if format == 'woff':
                        fonts[font] = url
    else:
        #print(rule.type)
        pass

#pprint(fonts)

for fn in fonts:
    print(fonts[fn])
    dfn = fn.replace(' ','_')
    print(''.join(["\tout=out/",dfn,'.woff.gz']))

'''
for rule in sheet:
    if rule.type == rule.FONT_FACE_RULE:
        # find property
        for property in rule.style:
            if property.name == 'font-family':
                font = property.value.strip('"')
                fonts[font] = {}
            if property.name == 'src':
                for pair in re.findall('url\((.*?)\) format\("(.*?)"\)', property.value):
                    #fonts[font][pair[1]] = pair[0]
                    if pair[1] == 'woff':
                        fonts[font] = pair[0]

uri可以读，但是cssfunction还是得regex，不如直接全部regex解决……

            if property.name == 'src':
                for i in range(0, property.propertyValue.length, 2):
                    url = property.propertyValue[i].absoluteUri
                    format = re.findall('"(.*?)"', property.propertyValue[i + 1].value)[0]
                    fonts[font][format] = url

wtff2otf:
https://github.com/hanikesn/woff2otf/blob/master/woff2otf.py
https://github.com/samboy/WOFF

https://github.com/google/woff2


otf2otc: (https://github.com/fonttools/fonttools/issues/17)
https://github.com/adobe-type-tools/afdko/blob/develop/python/afdko/otf2otc.py
https://github.com/googlei18n/nototools/blob/master/nototools/ttc_utils.py
'''
