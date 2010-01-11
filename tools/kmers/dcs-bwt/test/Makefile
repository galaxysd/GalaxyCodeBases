# Copyright 2007 Google Inc.
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
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

tests: stringsort_test difference_cover_test prefix_doubling_test
	./stringsort_test
	./difference_cover_test
	./prefix_doubling_test

stringsort_test: stringsort_test.cc

difference_cover_test: difference_cover_test.cc ../src/difference_cover.o

prefix_doubling_test: prefix_doubling_test.cc

clean:
	rm *~ stringsort_test difference_cover_test prefix_doubling_test
