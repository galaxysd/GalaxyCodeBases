from optparse import OptionParser
import cmarkgfm
from cmarkgfm.cmark import Options as cmarkgfmOptions
import os

parser=OptionParser()
parser.add_option("-i","--in_file", dest="infile", help="Input infile name", metavar="infile")
parser.add_option("-o","--out_file", dest="outfile", help="Input outfile name", metavar="outfile")
parser.add_option("-s","--style_file", dest="stfile", help="Input stfile name", metavar="stfile")
(options,args)=parser.parse_args()

Coptions = (
    cmarkgfmOptions.CMARK_OPT_UNSAFE
    | cmarkgfmOptions.CMARK_OPT_GITHUB_PRE_LANG
    | cmarkgfmOptions.CMARK_OPT_SMART
)

with open(options.infile,'r') as mdfile:
    html=cmarkgfm.github_flavored_markdown_to_html(mdfile.read(), Coptions)

htmlfile=open(options.outfile,"w")

with open(options.stfile,'r') as style:
    htmlfile.write(style.read())
htmlfile.write(html)
htmlfile.write("</body></html>")
htmlfile.close()
mdfile.close()
