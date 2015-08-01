# List

Buildcrx
========

From https://github.com/dangoakachan/buildcrx

This tool is to package chrome extension automatically, with below extra features:

* Exclude files or directories with the same pattern gramma like `grep -E`.
* Search manifest.json for update infomation, and generates the updates.xml automatically.

Package function is provided by [crxmake.sh](http://code.google.com/p/chromium/issues/attachmentText?id=15059&aid=-2305436989939443553&name=crxmake.sh).

Usage
=====

Script help:

    $ sh buildcrx.sh
    Usage: buildcrx.sh [OPTIONS]

    Options:
      -d     specify the extension directory, required
      -p     specify the pem file, optional
      -o     specify the output directory, optional
      -n     specify the build extension name, optional
      -e     specify the exclude pattern like grep, optional
      -h     show this help message

Example:

    $ sh buildcrx.sh -d copybetter -p copybetter.pem -o copybetter/update -e 'update|\.git'
