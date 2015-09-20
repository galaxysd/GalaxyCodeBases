#!/bin/bash
#
# This is a shell script to find the 'RealBin' and 'RealScript' for the
# given input command
#
# Change History:
#  1.0.1    - Basic release.
#



# Prints the absolute path of the given input directory.
# Usage: abs=`abs_path <dirName>`
function abs_path() {
    local dir=$1
    "cd" "$dir"
    if [ "$?" = "0" ]; then
        /bin/pwd
        "cd" -  &>/dev/null
    fi
}

# Finds the file from the 'PATH' environmental variable
# and then prints it's dir name.
# Usage: fromPATH <file>
function fromPATH() {
    local oldIfs=$IFS
    local file=$1
    local temp=""
    IFS=":"
    for dir in $PATH; do
        temp=$dir/$file
        if [ -r "$temp" ]; then
            temp=`/usr/bin/dirname "$temp"`
            abs_path "$temp"
            break
        fi
    done
    IFS=$oldIfs
}

# Equivalent to perl's FindBin::RealBin
# Usage: RealBin <file>
function RealBin() {
    local file=$1
    local temp=""
    if [ -r "$file" ]; then
        temp=`/usr/bin/dirname "$file"`
        abs_path "$temp"
    else
        fromPATH $file
    fi
}

# Equivalaent to perl's FindBin::RealScript
# Usage: RealScript <file>
function RealScript() {
    local file=$1
    local temp=`RealBin "$file"`
    temp=$temp/`/bin/basename "$file"`
    temp=`/usr/bin/readlink -f "$temp"`
    if [ "$temp" != "" ]; then
        /bin/basename "$temp"
    fi
}

# help and exit
function showHelp() {
    echo " FindBin.sh - Bash equivalent for Perl's FindBin.
USAGE:
  FindBin.sh [-h, -v] [-bin | -script] <file>
   -h        Print this help and exit.
   -v        Print version information and exit.
   -bin      Print the value of FindBin::RealBin.
   -script   Print the value of FindBin::RealScript.
   <file>    The file whose FindBin::RealBin & FindBin::RealScript
             are needed."
    exit 0
}

# version info
function showVersion() {
    echo "FindBin.sh - v1.0.1"
    exit 0
}


cmd="RealBin"
if [ "$1" = "-h" ]; then
    showHelp
elif [ "$1" = "-v" ]; then
    showVersion
elif [ "$1" = "-bin" ]; then
    cmd="RealBin"
    shift
elif [ "$1" = "-script" ]; then
    cmd="RealScript"
    shift
fi
"$cmd" "$1"
