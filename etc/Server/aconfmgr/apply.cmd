#!/bin/sh

mkdir files
ln -sf ../root-overlay files/
ln -sf ../calamares files/

#useradd -m galaxy
#usermod -aG wheel galaxy

#sudo ntpdate -u time.asia.apple.com
#sudo hwclock --systohc

# rsync -ahXUH 192.168.88.189:~/.config/aconfmgr/ . --delete -vrn
# rsync -ahXUH . 192.168.88.189:~/.config/aconfmgr/  -vn
