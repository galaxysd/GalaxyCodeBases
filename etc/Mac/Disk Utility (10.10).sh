#!/bin/sh -x ## or just ` curl -Ls http://git.io/vRozn | sh `.

## Downloads the Mac OS X 10.10 Recovery Partition update,
## Copy's over the 10.10 version of Disk Utility.app, then
## use git to apply a binary patch so it will run on 10.11+. 

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
