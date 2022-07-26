# set time
#CopyFile /etc/adjtime
if [[ $aconfmgr_action == "apply" ]]; then
	sudo /usr/bin/ntpdate -u time.asia.apple.com
	sudo /usr/bin/hwclock --systohc
fi

# set locale
CreateLink /etc/localtime /usr/share/zoneinfo/Asia/Hong_Kong
#CopyFile /etc/locale.gen
#CopyFile /etc/locale.conf

# dns
#CopyFile /etc/hostname
#CopyFile /etc/hosts

# pacman
#CopyFile /etc/pacman.conf
#CopyFile /etc/pacman.d/mirrorlist

# ssh
AddPackage openssh
AddPackage openssh-dinit

# sudo wheel
echo '%wheel ALL=(ALL) ALL' > "$(CreateFile /etc/sudoers.d/wheel)"

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

