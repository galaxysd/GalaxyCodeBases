#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import filecmp

dir1 = '''/Applications/Microsoft Word.app/'''
dir2 = '''/Applications/Microsoft Excel.app/'''


def get_common_files_helper(base_dir, dcmp):
    if os.path.islink(os.path.join(dcmp.right)):
        # print "skip dir link %s" % os.path.join(dcmp.left, base_dir)
        return []
    else:
        rs = []
        rs0 = filecmp.cmpfiles(dcmp.left, dcmp.right, dcmp.same_files, False)
        for name in rs0[0]:
            rs.append(os.path.join(base_dir, name))
        # for name in rs0[1]:
        #     print "%s is false same" % os.path.join(base_dir, name)

        for name, sub_dcmp in dcmp.subdirs.iteritems():
            path = os.path.join(base_dir, name)
            rs += get_common_files_helper(path, sub_dcmp)

        return rs


def get_common_files(dir1, dir2):
    return get_common_files_helper("", filecmp.dircmp(dir1, dir2))


def get_ino(path):
    return os.stat(path).st_ino

common_files = []
try:
    print "calculating common files..."
    common_files += get_common_files(dir1, dir2)

except OSError as ex:
    print "Error: %s" % ex.strerror
    print "Maybe you have installed Office in a non-default location."


n_common_files = len(common_files)
print "n_common_files: %s" % n_common_files
i = 0
for filename in common_files:
    path1 = os.path.join(dir1, filename)
    path2 = os.path.join(dir2, filename)
    # if get_ino(path1) == get_ino(path2):
    #     print "same inode for %s" % filename
    # else:
    #     print "diff inode for %s" % path1
    #     print "           and %s" % path2

    # os.rename won't create parent dir for me
    # os.renames will prune empty dir from src
    # Let's be lazy to skip backuping...
    os.remove(path2)
    os.link(path1, path2)

    i += 1
    if i % (n_common_files / 10) == 0:
        print '%s files linked...' % i

print "done"
