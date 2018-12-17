#!/bin/sh

FULLDATE=`date '+%Y%m%d%H%M%S'`
SSHCOMMENT="$(whoami)@$(hostname)-$(date '+%Y%m%d')"

ssh-keygen -b 4096 -t rsa -C ${SSHCOMMENT} -N '' -f id_rsa.4k.${FULLDATE}
ssh-keygen -b 3072 -t rsa -C ${SSHCOMMENT} -N '' -f id_rsa.3k.${FULLDATE}
ssh-keygen -t ecdsa -b 521 -C ${SSHCOMMENT} -N '' -f id_ecdsa.57.${FULLDATE}
ssh-keygen -t ed25519 -C ${SSHCOMMENT} -N '' -f id_ed25519.65.${FULLDATE}

ls -l id_*.${FULLDATE}*
echo id_*.${FULLDATE} | xargs -n1 ssh-keygen -lf
echo id_*.${FULLDATE}.pub | xargs -n1 ssh-keygen -lf

# https://www.keylength.com/
# https://wiki.archlinux.org/index.php/SSH_keys
# https://blog.josefsson.org/2016/11/03/why-i-dont-use-2048-or-4096-rsa-key-sizes/
# https://medium.com/@honglong/%E9%81%B8%E6%93%87-ssh-key-%E7%9A%84%E5%8A%A0%E5%AF%86%E6%BC%94%E7%AE%97%E6%B3%95-70ca45c94d8e
# echo 4096 | ./keysize-NIST.bc
