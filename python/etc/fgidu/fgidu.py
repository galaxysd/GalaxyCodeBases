#!/usr/bin/env python3

import os
# https://stackoverflow.com/a/62000501/159695
def du(path):
    if os.path.islink(path):
        return (os.lstat(path).st_size, 0)
    if os.path.isfile(path):
        st = os.lstat(path)
        return (st.st_size, st.st_blocks * 512)
    apparent_total_bytes = 0
    total_bytes = 0
    have = []
    for dirpath, dirnames, filenames in os.walk(path):
        apparent_total_bytes += os.lstat(dirpath).st_size
        total_bytes += os.lstat(dirpath).st_blocks * 512
        for f in filenames:
            fp = os.path.join(dirpath, f)
            if os.path.islink(fp):
                apparent_total_bytes += os.lstat(fp).st_size
                continue
            st = os.lstat(fp)
            if st.st_ino in have:
                continue  # skip hardlinks which were already counted
            have.append(st.st_ino)
            apparent_total_bytes += st.st_size
            total_bytes += st.st_blocks * 512
        for d in dirnames:
            dp = os.path.join(dirpath, d)
            if os.path.islink(dp):
                apparent_total_bytes += os.lstat(dp).st_size
    return (apparent_total_bytes, total_bytes)

def humanized_size(num, suffix='B', si=False, decps=1):
    # decps for decimal places
    if si:
        units = ['','K','M','G','T','P','E','Z']
        last_unit = 'Y'
        div = 1000.0
    else:
        units = ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']
        last_unit = 'Yi'
        div = 1024.0
    for unit in units:
        if abs(num) < div:
            #return "%3.1f%s%s" % (num, unit, suffix)
            return f"%3.{decps}f {unit}{suffix}" % (num)
        num /= div
    #return "%.1f%s%s" % (num, last_unit, suffix)
    return f"%.{decps}f {last_unit}{suffix}" % (num)

#s1,s2 = du("/Users/galaxy/t/t/t/t/Received")
s1,s2 = du("x")


print(s1)
print(s2)
print(humanized_size(s1,decps=2))
print(humanized_size(s2,si=True,decps=4))
