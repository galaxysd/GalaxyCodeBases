#!/usr/bin/env bash

FULLDATE=`date '+%Y%m%d%H%M%S'`
JUSTDATE=$(date '+%Y%m%d')
PLACES="local,bgi,vps"
# https://stackoverflow.com/a/45201229/159695
if [ "${BASH_VERSINFO:-0}" -ge 4 ]; then
	readarray -td',' aPLACES <<<"${PLACES},";unset aPLACES[-1]
else
	IFS=$',' read -ra aPLACES <<<"${PLACES}"
fi
#declare -p aPLACES; echo ${aPLACES[@]}
#echo $BASH_VERSION,$BASH_VERSINFO

for thePlace in ${aPLACES[@]}; do
	PREFIX="c${FULLDATE}.${thePlace}"
	echo ${PREFIX}
	mkdir ${PREFIX}
	ssh-keygen -q -b 4096 -t rsa -C " $(whoami)@${thePlace}.rsa4k.${JUSTDATE}" -N '' -f ${PREFIX}/id_rsa
	ssh-keygen -q -t ecdsa -b 521 -C " $(whoami)@${thePlace}.ecdsa521.${JUSTDATE}" -N '' -f ${PREFIX}/id_ecdsa
	ssh-keygen -q -t ed25519 -C " $(whoami)@${thePlace}.ed25519.${JUSTDATE}" -N '' -f ${PREFIX}/id_ed25519
	#cat ${PREFIX}/id_{ed25519,ecdsa,rsa}.pub >authorized_keys.${thePlace}.${FULLDATE}
done

rm -f authorized_keys.${FULLDATE}
for theType in rsa ecdsa ed25519; do
	rm -f authorized_keys.${theType}_${FULLDATE}
	for thePlace in ${aPLACES[@]}; do
		PREFIX="c${FULLDATE}.${thePlace}"
		cat ${PREFIX}/id_${theType}.pub >>authorized_keys.${theType}_${FULLDATE}
	done
	cat authorized_keys.${theType}_${FULLDATE} >>authorized_keys.${FULLDATE}
done
