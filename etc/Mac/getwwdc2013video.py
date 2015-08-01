#!/usr/bin/env python
# coding: utf-8
__author__ = 'Snake'


import urllib2, cookielib
from bs4 import BeautifulSoup as bs

base_url = 'https://developer.apple.com/wwdc/videos/'
base_header = {'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; rv:21.0) Gecko/20100101 Firefox/21.0'}


#获取网站内容
def get_web(url, cookie, choose):
    download = []
    content = get_url_result(url, cookie)
    soup = bs(str(content))
    for div in soup.findAll('p', {'class': 'download'}):
        if choose != 1 and choose != 3 and choose != 4:
            hd =  div('a')[0]['href']
            download.append(hd)
        if choose != 2 and choose != 3 and choose != 5:
            sd = div('a')[1]['href']
            download.append(sd)
        try:
            if choose != 1 and choose != 2:
                pdf = div('a')[2]['href']
                download.append(pdf)
        except:
            print "none pdf"
    return download


#根据url和cookie获取网页内容
#如果传入cookie不为空那么将通过cookie的方式获取url信息
#如果传入的cookie为空那么就已普通的方式过去url信息
def get_url_result(url, cookie=None):
    headers = base_header
    headers['Referer'] = url
    if cookie:
        headers['Cookie'] = cookie
        opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cookielib.CookieJar()))
        urllib2.install_opener(opener)
    req = urllib2.Request(url=url, headers=headers)
    try:
        content = urllib2.urlopen(req, timeout=20).read()
        return content
    except  Exception as e:
        print e.message + '\n'



def main():
    cookie = raw_input("cookie: ")
    print "please choose: 1-SD,2-HD,3-PDF,4-SD+PDF,5-HD+PDF,0-ALL"
    choose = int(raw_input("choose: "))
    download = get_web(base_url, cookie, choose)
    for d in download:
        print d

if __name__ == '__main__':
    main()
