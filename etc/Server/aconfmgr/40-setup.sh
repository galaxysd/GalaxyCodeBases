# set time
IgnorePath '/etc/adjtime'
if [[ $aconfmgr_action == "apply" ]]; then
	sudo /usr/bin/ntpdate -u time.asia.apple.com || true
	sudo /usr/bin/hwclock --systohc
fi
CreateLink /etc/localtime /usr/share/zoneinfo/Asia/Hong_Kong

# dns
#CopyFile /etc/hostname
#CopyFile /etc/hosts

# pacman
#CopyFile /etc/pacman.conf
#CopyFile /etc/pacman.d/mirrorlist

# sudo wheel
echo '%wheel ALL=(ALL) ALL' > "$(CreateFile /etc/sudoers.d/wheel 440)"

# Add user
if [[ $aconfmgr_action == "apply" ]]; then
	if id -u "galaxy" >/dev/null 2>&1; then
		echo "[!]User Galaxy exists."
	else
		sudo useradd -m galaxy
		sudo usermod -aG wheel galaxy
	fi
fi

AddPackage expac
#AddPackage lesspipe
AddPackage noto-fonts-cjk
AddPackage noto-fonts-emoji
AddPackage noto-fonts-extra
AddPackage git
AddPackage gnuplot
AddPackage python-pygments
AddPackage bash-completion
AddPackage maliit-keyboard
AddPackage python-markdown
AddPackage poppler-data
#AddPackage gocryptfs
AddPackage unarchiver
AddPackage htop
#AddPackage iotop

AddPackage --foreign aconfmgr-git # A configuration manager for Arch Linux
AddPackage --foreign yay # Yet another yogurt. Pacman wrapper and AUR helper written in go.

AddPackage vim
AddPackage lsof
#AddPackage strace
AddPackage ed

AddPackage linux-zen
AddPackage linux-zen-headers
AddPackage go

