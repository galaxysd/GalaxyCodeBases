#!/bin/bash

function pause(){
	tty_state=$(stty -g)
	stty -icanon min 0 time 0
	while read dummy; do : ; done
	stty "$tty_state"
	read -n1 -rsp "$*" < /dev/tty
	echo
}

INSTPATH=/Users/Shared/Fonts/Noto/Noto-unhinted

mkdir -p $INSTPATH

if [ x"$1" = x"" ]; then
	echo "Install downloaded zip files ..."
	if [[ ! -f 'NotoSansCJK.ttc.zip' ]]; then
		echo "[NotoSansCJK.ttc.zip] not found !"
		exit 1
	fi
	if [[ ! -f 'Noto-unhinted.zip' ]]; then
		echo "[Noto-unhinted.zip] not found !"
		exit 1
	fi
else
	echo "Downloading Google files ..."
	mv NotoSansCJK.ttc.zip NotoSansCJK.ttc.zip.last
	mv Noto-unhinted.zip Noto-unhinted.zip.last

	# http://www.google.com/get/noto/help/cjk/
	wget https://noto-website-2.storage.googleapis.com/pkgs/NotoSansCJK.ttc.zip
	wget https://noto-website-2.storage.googleapis.com/pkgs/Noto-unhinted.zip
fi

ls -l NotoSansCJK.* Noto-unhinted.*
shasum NotoSansCJK.* Noto-unhinted.*

pause 'Press any key to continue or Ctrl+C to exit ... '

echo "Unzip files ..."
unzip -qod /Users/Shared/Fonts/Noto NotoSansCJK.ttc.zip
unzip -qod ${INSTPATH} Noto-unhinted.zip

ls -1 ${INSTPATH}/NotoSans*CJK* |cat -n

pause $'There shoud be 36 files.\nPress any key to continue ... '

rm ${INSTPATH}/NotoSans*CJK*
find /Users/Shared/Fonts/Noto -name '*.txt' -delete
echo "Cleaning done. You can remove NotoSansCJK.* Noto-unhinted.* now."

open /Users/Shared/Fonts/Noto
