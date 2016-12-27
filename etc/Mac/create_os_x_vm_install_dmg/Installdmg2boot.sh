#!/bin/bash
#
# This executable converts a Mavericks .app (which allows to upgrade a machine
# from Mac OS 10.6.7+ to Mavericks) into a Mavericks .dmg (which allows to
# install Mavericks from scratch on a machine).
#
# It has been tested with "Install OS X 10.9 Developer Preview.app" (build
# 13A476u).
#
# http://www.insanelymac.com/forum/topic/290949-how-to-install-os-x-10x-snow-leopard-to-el-capitan-in-vmware-workstation-1011-workstation-proplayer-12-player-67-esxi-56/page-2#entry1956505

set -x
set -e


# The first argument is the path to the .app bundle (the input of the
# executable).
inputApp="$1"
# The second argument is the path to the .dmg file (the output of the
# executable), which must end with ".dmg".
outputDmg="$2"
[ "${outputDmg: -4}" = .dmg ]


#
# The problem: /System/Installation/Packages inside /BaseSystem.dmg inside
# "$inputApp"/Contents/SharedSupport/InstallESD.dmg is a dangling symlink,
# which prevents installing Mavericks from scratch.
# The solution: Replace the symlink with the /Packages directory inside
# "$inputApp"/Contents/SharedSupport/InstallESD.dmg.
#


tmpDir=`mktemp -d -t 'Create Mavericks Installer'`
installMnt="$tmpDir"/install
installPkg="$installMnt"/Packages
outputMnt="$tmpDir"/output
outputPkg="$outputMnt"/System/Installation/Packages


cleanup() {
   if [ -d "$outputMnt" ]; then
      hdiutil detach "$outputMnt"
   fi


   if [ -d "$installMnt" ]; then
      hdiutil detach "$installMnt"
   fi


   rmdir -- "$tmpDir"
}


# Cleanup on failure.
trap cleanup ERR


# Mount InstallESD.dmg so we can access /BaseSystem.dmg and /Packages inside.
hdiutil attach "$inputApp"/Contents/SharedSupport/InstallESD.dmg \
   -mountpoint "$installMnt" -nobrowse


# Create "$outputDmg", a read/write copy of the read-only BaseSystem.dmg.
hdiutil convert "$installMnt"/BaseSystem.dmg -format UDRW -o "$outputDmg"


# Enlarge "$outputDmg" to accommodate for our modifications. The UDRW image
# format is not sparse, so we must precisely compute the new size.
curSectors=`hdiutil resize "$outputDmg" -limits | tail -1 | awk '{ print $2 }'`
extraSectors=`BLOCKSIZE=512 du -s -- "$installPkg" | awk '{ print $1 }'`
hdiutil resize "$outputDmg" -sectors $((curSectors + extraSectors))


# Mount "$outputDmg".
hdiutil attach "$outputDmg" -mountpoint "$outputMnt" -nobrowse


# Modify "$outputDmg".
rm -- "$outputPkg"
cp -r -- "$installPkg" "$outputPkg"


# Cleanup on success.
trap ERR; cleanup


ls -alh -- "$outputDmg"
