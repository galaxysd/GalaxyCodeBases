f="$(GetPackageOriginalFile grub /etc/default/grub)"
sed -E -i 's/^#?(GRUB_TIMEOUT)=.*/\1=3/g' "$f"
sed -E -i 's/^#?(GRUB_CMDLINE_LINUX_DEFAULT)=.*/\1="quiet splash"/g' "$f"
sed -E -i 's/^#?(GRUB_CMDLINE_LINUX)=.*/\1="net.ifnames=0"/g' "$f"
sed -E -i 's/^#?(GRUB_GFXMODE)=.*/\1="1280x800,1024x768,800x600"/g' "$f"
sed -E -i 's|^#?(GRUB_THEME)=.*|\1="/usr/share/grub/themes/artix/theme.txt"|g' "$f"

cat >> "$(GetPackageOriginalFile pam /etc/environment)" <<EOF
QT_STYLE_OVERRIDE=gtk2
QT_QPA_PLATFORMTHEME=gtk
QTWEBENGINE_CHROMIUM_FLAGS="-blink-settings=darkModeEnabled=true -enable-features=OverlayScrollbar,OverlayScrollbarFlashAfterAnyScrollUpdate,OverlayScrollbarFlashWhenMouseEnter"
EOF

cat >> "$(GetPackageOriginalFile filesystem /etc/hosts)" <<EOF
127.0.0.1	localhost
::1		localhost ip6-localhost ip6-loopback
ff02::1		ip6-allnodes
ff02::2		ip6-allrouters
EOF

CopyFileTo "root-overlay/etc/sddm.conf" "/etc/sddm.conf"
CopyFileTo "root-overlay/usr/share/gtk-2.0/gtkrc" "/usr/share/gtk-2.0/gtkrc"

cat > "$(CreateFile /etc/udev/rules.d/jms580583trim.rules)" <<EOF
ACTION=="add|change", ATTRS{idVendor}=="152d", ATTRS{idProduct}=="a583", SUBSYSTEM=="scsi_disk", ATTR{provisioning_mode}="unmap", ATTR{manage_start_stop}="1"
EOF

# Specify locales
f="$(GetPackageOriginalFile glibc /etc/locale.gen)"
sed -i 's/^#\(en_US.UTF-8\)/\1/g' "$f"
sed -i 's/^#\(en_HK.UTF-8\)/\1/g' "$f"
sed -i 's/^#\(zh_CN.UTF-8\)/\1/g' "$f"
sed -i 's/^#\(zh_TW.UTF-8\)/\1/g' "$f"
sed -i 's/^#\(ja_JP.UTF-8\)/\1/g' "$f"

if [[ $aconfmgr_action == "apply" ]]; then
	sudo /usr/bin/locale-gen
fi

