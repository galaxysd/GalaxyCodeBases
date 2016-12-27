#!/bin/bash
# http://www.insanelymac.com/forum/topic/293481-how-to-create-a-bootable-iso-from-the-mavericks-installesddmg/
# 测试后发现无法启动……

ESD=$1
TMP=$2

if [ -z "$ESD" ] || [ -z "$TMP" ]; then
	echo usage: "'$0' /path/to/esd /path/to/tmpdir"
	exit 1
fi
if ! [ -e "$ESD" ]; then
	echo "file '$ESD' does not exist"
	exit 1
fi
if ! [ -e "$TMP" ]; then
	echo "dir '$TMP' does not exist"
	exit 1
fi

MPAPP=/Volumes/install_app
MPIMG=/Volumes/install_img
IMGSPARSE=$TMP/install.sparseimage
IMGDVD=$TMP/install.cdr

detach_all() {
  if [ -d "$MPAPP" ]; then hdiutil detach "$MPAPP"; fi
  if [ -d "$MPIMG" ]; then hdiutil detach "$MPIMG"; fi
}
exit_all() {
  echo +++ Command returned with error, aborting ...
  exit 2
}

trap detach_all EXIT
trap exit_all ERR

echo +++ Trying to unmount anything from previous run
detach_all

echo +++ Mount the installer image
hdiutil attach "$ESD" -noverify -nobrowse -readonly -mountpoint "$MPAPP"

echo +++ Convert the boot image to a sparse bundle
rm -f "$IMGSPARSE"
hdiutil convert "$MPAPP"/BaseSystem.dmg -format UDSP -o "$IMGSPARSE"

echo +++ Increase the sparse bundle capacity to accommodate the packages
hdiutil resize -size 9g "$IMGSPARSE"

echo +++ Mount the sparse bundle for package addition
hdiutil attach "$IMGSPARSE" -noverify -nobrowse -readwrite -mountpoint "$MPIMG"

echo +++ Remove Package link and replace with actual files
rm -f "$MPIMG"/System/Installation/Packages
cp -rp "$MPAPP"/Packages "$MPIMG"/System/Installation/

echo +++ Unmount the installer image
hdiutil detach "$MPAPP"

echo +++ Unmount the sparse bundle
hdiutil detach "$MPIMG"

echo +++ Resize the partition in the sparse bundle to remove any free space
hdiutil resize -sectors min "$IMGSPARSE"

echo +++ Convert the sparse bundle to ISO/CD master
rm -f "$IMGDVD"
hdiutil convert "$IMGSPARSE" -format UDTO -o "$IMGDVD"

echo +++ Remove the sparse bundle
rm "$IMGSPARSE"


echo "Done"
echo "Find your DVD at '$IMGDVD'"
