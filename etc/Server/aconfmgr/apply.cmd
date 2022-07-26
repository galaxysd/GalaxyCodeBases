#!/bin/sh

mkdir files
ln -sf ../root-overlay files/

useradd -m galaxy
usermod -aG wheel galaxy

sudo ntpdate -u time.asia.apple.com
sudo hwclock --systohc

