#!/bin/sh

git remote add g git@gitorious.org:galaxycodebases/main.git
git remote add c git@coding.net:galaxy/GalaxyCodeBase.git
git remote -v
git remote update

cd released/pIRS
git remote add github git@github.com:galaxy001/pirs.git
git remote add google https://code.google.com/p/pirs/
git remote -v
git remote update

git config --global alias.xpush '!f(){ for i in `git remote`; do git push -v $i; done; };f'
git config --global alias.xpull '!f(){ for i in `git remote`; do git pull -v $i; done; };f'

