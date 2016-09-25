#!/bin/sh -x ## or just ` curl -Ls http://git.io/vRozn | sh `.

## Downloads the Mac OS X 10.10 Recovery Partition update,
## Copy's over the 10.10 version of Disk Utility.app, then
## use git to apply a binary patch so it will run on 10.11+. 

# https://gist.github.com/geoff-codes/b96bc9c5bca538aa0819
# https://justus.berlin/2015/10/restore-old-disk-utility-in-os-x-el-capitan/

# http://www.insanelymac.com/forum/files/file/480-disk-utilitypatched/
# Disk Utility 10.10.5
# defaults write com.apple.DiskUtility DUDebugMenuEnabled 1
# defaults write com.apple.DiskUtility advanced-image-options 1

# http://www.insanelymac.com/forum/files/file/481-diskutilitypatch/
# systemVersionCheck
# perl -pi -e 's|\xD5\x84\xC0\x0F\x85\x44\x01\x00\x00\x48\x8B\x05\x10\x5B\x06\x00|\xD5\x84\xC0\xE9\x45\x01\x00\x00\x00\x48\x8B\x05\x10\x5B\x06\x00|g'
#
# enableDebugMenu 
# debugMode can also be achieved by drag and dropping the Disk Utility binary on a terminal window and adding --debugMenu at the end.
# perl -pi -e 's|\x3C\x33\x74\x48\x48\x85|\x3C\x33\x66\x90\x48\x85|g'
# perl -pi -e 's|\x85\xFF\x75\x58\x48\x8B|\x85\xFF\x66\x90\x48\x8B|g'
# 
# expertMode
# perl -pi -e 's|\xFF\x00\x75\x25\x4C\x89\xEF|\xFF\x00\xEB\x25\x4C\x89\xEF|g' [Drag & Drop Disk Utility binary]
# 
# md5 for original 10.10.5 DiskUtility v13 (606)	f4de87e306c1bc6acdef6c6166e574ca
# md5 for systemVersionCheck binary (only)	3890ef83d60c9e1838be6cf2eee35d18
# md5 combo systemVersioncheck & expertMode & enableDebugMenu	28e835bf915205c922b892cf01b5b19f

# After you have downloaded and unzipped the file, you need to let root:wheel take ownership of the app ie
# sudo chown -R 0:0 ~/Downloads/Disk\ Utility.app

cd /tmp

rm -rf DU1010
mkdir DU1010
cd DU1010

# See: https://github.com/wdas/reposado
curl -ORL http://swcdn.apple.com/content/downloads/21/09/031-20634/8d84o1ky5gn2agnf5kiz9eed134n7y3q4c/RecoveryHDUpdate.pkg

xar -xf RecoveryHDUpdate.pkg

hdiutil attach -nobrowse RecoveryHDMeta.dmg
hdiutil attach -nobrowse "/Volumes/Recovery HD Update/BaseSystem.dmg"

cp -Rp "/Volumes/OS X Base System/Applications/Utilities/Disk Utility.app" .

hdiutil detach "/Volumes/OS X Base System"
hdiutil detach "/Volumes/Recovery HD Update"

rm * 2>/dev/null

git apply --no-index <<-PATCH
diff --git a/Disk Utility.app/Contents/MacOS/Disk Utility b/Disk Utility.app/Contents/MacOS/Disk Utility
index 20f16f995699699ac5841a455c261a2e45e6a4e5..ce4994ac078a99803ba6b5f3467fae425b5bf643 100755
GIT binary patch
delta 44
scmZoTC);pNcEgiI#+K%%iS17l8G)Dyh?#+y1&CRJm~H#hM0ROb0OUIqApigX

delta 44
scmZoTC);pNcEgiI#@6PiiS17l8G)Dyh?#+y1&CRJm~H#hM0ROb0OVa1A^-pY

PATCH

mv "Disk Utility.app" "Disk Utility (10.10).app"

open .
