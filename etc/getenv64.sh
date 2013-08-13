#!/bin/sh
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#
# Copyright (c) 2006 Guillaume Chazarain <guichaz@gmail.com>
#
# History:
#  20061027: Initial release
#  20110314: Compile with -fPIC

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 COMMAND [ARGS...]" >&2
    exit 1
fi

TEMPDIR=$(mktemp -d)
[ -z "$TEMPDIR" ] && exit 1
trap "rm -fr $TEMPDIR" EXIT

cat > "$TEMPDIR/getenv.c" <<'EOF' || exit 1
#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>

char *getenv (const char *name)
{
	static char *(*real_getenv)(const char *name) = NULL;
	char *ret;

	if (real_getenv == NULL) {
		void *handle = dlopen("/lib64/libc.so.6", RTLD_LAZY);
		if (handle == NULL) {
			fprintf(stderr, "dlopen: %s\n", dlerror());
			exit(1);
		}

		real_getenv = dlsym(handle, "getenv");
		if (real_getenv == NULL) {
			fprintf(stderr, "dlsym: %s\n", dlerror());
			exit(1);
		}
	}

	ret = real_getenv(name);
	if (ret == NULL)
		printf("getenv(\"%s\") = NULL\n", name);
	else
		printf("getenv(\"%s\") = \"%s\"\n", name, ret);

	return ret;
}
EOF

gcc -fPIC -O0 -ldl -shared "$TEMPDIR/getenv.c" -o "$TEMPDIR/libgetenv.so" || exit 1
LD_PRELOAD="$TEMPDIR/libgetenv.so:$LD_PRELOAD" "$@"
