#!/bin/sh

# git clone git@github.com:galaxysd/GalaxyCodeBases.git toGit

git remote add gitlab git@gitlab.com:galaxycodebases/main.git
git remote add coding git@git.coding.net:galaxy/GalaxyCodeBase.git
git remote add orcz ssh://galaxy@repo.or.cz/GalaxyCodeBases.git
# git remote add github git@github.com:galaxysd/GalaxyCodeBases.git
# git push -u github master
git remote -v
git remote update

git submodule init
git submodule update

git submodule foreach git branch -v
git submodule foreach git checkout master
git submodule foreach git pull
git submodule foreach git branch -v

cd released/pIRS
git remote add sf ssh://galaxy001@git.code.sf.net/p/pirsim/code
# git remote add github git@github.com:galaxy001/pirs.git
# git remote add google https://code.google.com/p/pirs/
git remote -v
git remote update

git config --global alias.xpush '!f(){ for i in `git remote`; do git push -v $i; done; };f'
git config --global alias.xpull '!f(){ for i in `git remote`; do git pull -v $i $1; done; };f'

git config --global core.editor "vim"

# git submodule add git@bitbucket.org:galaxy001/mydoc.git mydoc
#git submodule set-url mydoc git@bitbucket.org:galaxysd/mydoc.git
