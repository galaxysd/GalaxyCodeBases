

# Tuesday, July 26, 2022 PM04:05:33 HKT - Unknown packages


AddPackage wireless-regdb # Central Regulatory Domain Database


# Tuesday, July 26, 2022 PM04:05:33 HKT - Missing packages


RemovePackage crda


# Tuesday, July 26, 2022 PM04:05:33 HKT - Unknown foreign packages


AddPackage --foreign aconfmgr-git # A configuration manager for Arch Linux


# Tuesday, July 26, 2022 PM04:05:33 HKT - New / changed files


CreateDir /var/lib/libvirt/images
CreateDir /var/lib/rpcbind 700 rpc rpc
CreateDir /var/lib/sddm '' sddm sddm
CreateDir /var/lib/tpm2-tss/system/keystore 2775 tss tss


# Tuesday, July 26, 2022 PM04:05:33 HKT - New file properties


SetFileProperty /etc/sudoers.d/wheel mode 440
SetFileProperty /usr/bin/newgidmap mode 755
SetFileProperty /usr/bin/newuidmap mode 755
SetFileProperty /usr/lib/utempter/utempter group utmp
SetFileProperty /usr/lib/utempter/utempter mode 2755
