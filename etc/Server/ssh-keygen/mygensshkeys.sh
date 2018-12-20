#!/usr/bin/env bash

# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash/14203146#14203146
function show_help {
	echo "Usage: $0 [-hcv]"
	echo "Options: -h   display this help"
	echo "         -c   clean previous key files"
	echo "         -v   be verbose"
}
function do_clean {
	rm -rv c* authorized_keys.*
}

OPTIND=1
output_file=""
verbose=0

while getopts "h?cvf:" opt; do
	case "$opt" in
	h|\?)
		show_help
		exit 0
		;;
	c)
		do_clean
		exit 0
		;;
	v)  verbose=1
		;;
	f)  output_file=$OPTARG
		;;
    esac
done

shift $((OPTIND-1))
[ "${1:-}" = "--" ] && shift
#echo "verbose=$verbose, output_file='$output_file', Leftovers: $@"

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

rm -f authorized_keys.${JUSTDATE}
for theType in rsa ecdsa ed25519; do
	rm -f authorized_keys.${theType}_${JUSTDATE}
	for thePlace in ${aPLACES[@]}; do
		PREFIX="c${FULLDATE}.${thePlace}"
		cat ${PREFIX}/id_${theType}.pub >>authorized_keys.${theType}_${JUSTDATE}
	done
	cat authorized_keys.${theType}_${JUSTDATE} >>authorized_keys.${JUSTDATE}
done
