#!/bin/bash

source_disk=/dev/disk0s4
host='huxs@lab.luo-lab.org'
host_image=/share/users/huxs/t/disk0s4.img

bytes_transferred() {
  if ! ssh -p2211 "$host" "test -e '$host_image'"; then
    echo 0
    return
  fi
  ssh -p2211 "$host" "stat -c '%s' '$host_image'"
}
bytes_total() {
  #echo $(( $(blockdev --getsz $source_disk) * 512 ))
  echo `diskutil info $source_disk|grep 'Total Size'|awk '{print $5}'|sed 's/(//'`
}

echo "$(bytes_total) to $(bytes_transferred)"

while (( $(bytes_transferred) < $(bytes_total) )); do
  ( dd bs=1 skip=$(bytes_transferred) count=0 2>/dev/null && cat ) < $source_disk |lz4 -9 - - | ssh -p2211 -T -c arcfour -o Compression=no -x "$host" "lz4cat - >> '$host_image'"
done

