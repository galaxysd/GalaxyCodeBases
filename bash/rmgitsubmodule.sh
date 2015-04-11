#!/bin/sh

# http://stackoverflow.com/questions/1260748/how-do-i-remove-a-git-submodule

if [ x"$1" = x"" ]; then
	echo "Usage: $0 <git submodule name>"
	exit 1
else
	echo "Removing submodule [$1] ..."
fi

git submodule deinit $1    
git rm $1

# Note: [asubmodule] (no trailing slash)
# or, if you want to leave it in your working tree
#git rm --cached $1

echo Cleanup [.git/modules/$1]
rm -rf .git/modules/$1
