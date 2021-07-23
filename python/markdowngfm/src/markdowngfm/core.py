# Insert your code here. 
from weasyprint import HTML
import datetime
import pkg_resources
import os

import cmarkgfm
from cmarkgfm.cmark import Options as cmarkgfmOptions
_GFMoptions_ = (
    cmarkgfmOptions.CMARK_OPT_GITHUB_PRE_LANG
    | cmarkgfmOptions.CMARK_OPT_SMART
    | cmarkgfmOptions.CMARK_OPT_UNSAFE
)   # CMARK_OPT_UNSAFE is required for inline image/svg+xml

import pprint
pp = pprint.PrettyPrinter(indent=4)

def mainfunc(nfofile,msfile,hlafile,medfiles,outpath):
    """Working script."""
    print('Main')
    print(nfofile)
    print(msfile)
    print(hlafile)
    print(medfiles)
    print(outpath)

def getdat(resname):
    stream = pkg_resources.resource_stream(__name__, resname)
    read_data = stream.read().decode("utf-8", "ignore")
    return read_data

def _write2txt(mdstr,fm,fh):
    html = cmarkgfm.github_flavored_markdown_to_html(mdstr, _GFMoptions_)
    print(mdstr, file=fm)
    print(html, file=fh)
