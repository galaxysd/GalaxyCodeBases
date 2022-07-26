# core packages
AddPackage base
AddPackageGroup base-devel

# dinit
AddPackage dinit
AddPackage elogind-dinit

# kernel
AddPackage linux-zen
AddPackage linux-firmware
AddPackage linux-zen-headers

# filysystem
#CopyFile /etc/fstab

# system-boot
AddPackage efibootmgr
AddPackage grub
AddPackage os-prober
#CopyFile /etc/machine-id 444
#CreateLink /var/lib/dbus/machine-id /etc/machine-id

# initramfs
#CopyFile /etc/mkinitcpio.d/linux-zen.preset
#CopyFile /etc/arch-release
#CopyFile /etc/artix-release
CreateLink /etc/os-release ../usr/lib/os-release

