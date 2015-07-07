#!/bin/sh
#
# http://modmyi.com/forums/general/713564-how-repack-dpkgs.html#post7205248
#
# redeb.sh - Command line utility for rewrapping installed Debian-packages
# Copyright (C) 2015 Ikem Krueger <ikem.krueger@gmail.com>
# Copyright (C) 2003 Tommi Saviranta <tsaviran@cs.helsinki.fi>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
#
# Version: redeb.sh v1.0.0  15-January-2015  ikem.krueger@gmail.com

BASENAME="${0##*/}"
PWD_OLD="$PWD"
PACKAGE="$1"
TEMP="/tmp/$BASENAME-$$/$PACKAGE"

if [ $# -ne 1 ]
then
	echo "Usage: $BASENAME <PACKAGE>"

	exit 0
fi

# See if such package is installed.
if [ ! -f "/var/lib/dpkg/info/$PACKAGE.list" ]
then
	echo "No such package installed."

	exit 1
fi

# Create temporary directory.
mkdir -p "$TEMP"
cd "$TEMP"

# Get file list for this package.
for FILE in $(cat "/var/lib/dpkg/info/$PACKAGE.list")
do
	test -f "$FILE" && echo "$FILE">>filelist.tmp
done

# Copy files here.
tar -cf - $(cat filelist.tmp)|tar -xf -

# Create package-info.
mkdir DEBIAN

cp /var/lib/dpkg/info/$PACKAGE.* DEBIAN

rm "DEBIAN/$PACKAGE.list"

TOTAL=$(wc -l /var/lib/dpkg/available|cut -d ' ' -f 1)
BEGIN=$(grep -nx "Package: $PACKAGE" /var/lib/dpkg/available|cut -d ':' -f 1)

tail -n $(expr $TOTAL - $BEGIN + 1) /var/lib/dpkg/available>control.tmp

LENGTH=$(grep -nx "" control.tmp|cut -d ':' -f 1|head -n 1)

head -n $LENGTH control.tmp>DEBIAN/control

# Get package info
VERSION=$(grep -E "^Version:" DEBIAN/control|cut -d ' ' -f 2|cut -d ':' -f 2)
ARCH=$(grep -E "^Architecture:" DEBIAN/control|cut -d ' ' -f 2)

# Clean up.
rm -f control.tmp filelist.tmp

# Create package.
cd ..
dpkg-deb -b "$TEMP" "$PWD_OLD/${PACKAGE}_${VERSION}_${ARCH}.deb"

# Clean up the rest.
rm -rf "$TEMP"
