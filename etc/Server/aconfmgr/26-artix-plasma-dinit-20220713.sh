AddPackage cups-dinit # dinit service scripts for cups
CreateLink /etc/dinit.d/boot.d/cupsd /etc/dinit.d/cupsd

AddPackage networkmanager-dinit # dinit service scripts for networkmanager
CreateLink /etc/dinit.d/boot.d/NetworkManager /etc/dinit.d/NetworkManager

AddPackage sddm-dinit # dinit service scripts for sddm
CreateLink /etc/dinit.d/boot.d/sddm /etc/dinit.d/sddm

AddPackage syslog-ng-dinit # dinit service scripts for syslog-ng
CreateLink /etc/dinit.d/boot.d/syslog-ng /etc/dinit.d/syslog-ng

AddPackage wireless-regdb # Central Regulatory Domain Database

AddPackage openssh
AddPackage openssh-dinit
CreateLink /etc/dinit.d/boot.d/sshd ../sshd

CreateLink /etc/dinit.d/boot.d/acpid /etc/dinit.d/acpid
CreateLink /etc/dinit.d/boot.d/bluetoothd /etc/dinit.d/bluetoothd
CreateLink /etc/dinit.d/boot.d/cronie /etc/dinit.d/cronie

CopyFileTo "calamares/etc/default/locale" "/etc/default/locale"
CopyFileTo "calamares/etc/X11/xorg.conf.d/00-keyboard.conf" "/etc/X11/xorg.conf.d/00-keyboard.conf"
CopyFileTo "calamares/etc/default/keyboard" "/etc/default/keyboard"
#CopyFileTo "calamares/etc/hostname" "/etc/hostname"
CopyFileTo "calamares/etc/local.d/branding.start" "/etc/local.d/branding.start"
CopyFileTo "calamares/etc/locale.conf" "/etc/locale.conf"
CopyFileTo "calamares/etc/timezone" "/etc/timezone"
CopyFileTo "calamares/etc/vconsole.conf" "/etc/vconsole.conf"
CopyFileTo "calamares/etc/crypttab" "/etc/crypttab"

# from 99-unsorted.sh
CreateDir /var/lib/rpcbind 700 rpc rpc
CreateDir /var/lib/tpm2-tss/system/keystore 2775 tss tss
SetFileProperty /etc/cups/classes.conf mode 600
SetFileProperty /etc/cups/printers.conf mode 600
SetFileProperty /usr/bin/newgidmap mode 755
SetFileProperty /usr/bin/newuidmap mode 755
SetFileProperty /usr/lib/utempter/utempter group utmp
SetFileProperty /usr/lib/utempter/utempter mode 2755
