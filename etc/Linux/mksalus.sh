#!/bin/sh

mkdir salusSpatialOmics
cd salusSpatialOmics
git clone --bare --recursive --filter=blob:none --also-filter-submodules github-BGI:salusbio/SpatialOmics.git

cd SpatialOmics.git/
git config --add remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*"
git fetch

git config --local user.name 'HU Xuesong' 
git config --local user.email '87519979+huxs001@users.noreply.github.com'

git for-each-ref 'refs/heads/*' --format '%(refname:short)' --sort=-creatordate --shell | head -n5 | xargs -I{} echo "git worktree add ../{} {}; git branch -u origin/{} {}" | bash

git branch -av

# https://stackoverflow.com/a/74876298/159695 and https://stackoverflow.com/a/62524752/159695 
# for `git for-each-ref --format='%(refname:short)' refs/heads | xargs git branch -d`, there will be not 'refs/heads/*' but only 'refs/remotes/origin/*', and its short refname is like `origin/main`, which is not good.
# However, `git worktree add ../main main` works.
# https://stackoverflow.com/questions/14639206/how-can-i-pass-all-arguments-with-xargs-in-middle-of-command-in-linux/35612138#35612138
