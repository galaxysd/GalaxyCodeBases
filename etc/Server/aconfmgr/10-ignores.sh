# kernel
IgnorePath '/boot/vmlinuz-linux-zen'
IgnorePath '/boot/initramfs-linux-zen.img'
IgnorePath '/boot/initramfs-linux-zen-fallback.img'
IgnorePath '/etc/mkinitcpio.d/linux-zen.preset'
IgnorePath '/usr/lib/modules/*-zen*'

# UEFI
IgnorePath '/boot/efi/*'

# pacman
IgnorePath '/etc/pacman.d/gnupg'
IgnorePath '/var/lib/pacman/local/*' # package metadata
IgnorePath '/var/lib/pacman/sync/*.db' # repos
IgnorePath '/var/lib/pacman/sync/*.db.sig' # repo sigs

# certificates
IgnorePath '/etc/ca-certificates/extracted'
IgnorePath '/etc/ssl/certs'

# ld
IgnorePath '/etc/ld.so.cache'

# info
IgnorePath '/usr/share/info/dir'

# udev
IgnorePath '/etc/udev/hwdb.bin'
IgnorePath '/usr/lib/udev/hwdb.bin'

# system logs
IgnorePath '/var/log/*'

# sudo lectured
IgnorePath '/var/db/sudo/lectured/*'

# gtk
IgnorePath '/usr/lib/gdk-pixbuf-2.0/2.10.0/loaders.cache'
IgnorePath '/usr/lib/gtk-2.0/2.10.0/immodules.cache'
IgnorePath '/usr/lib/gtk-3.0/3.0.0/immodules.cache'
IgnorePath '/usr/lib32/gtk-2.0/2.10.0/immodules.cache'

# gnome
IgnorePath '/usr/lib/gio/modules/giomodule.cache'

# mime
IgnorePath '/usr/share/applications/mimeinfo.cache'
IgnorePath '/usr/share/mime/*.xml' # Localizations
IgnorePath '/usr/share/mime/XMLnamespaces'
IgnorePath '/usr/share/mime/aliases' # MIME aliases
IgnorePath '/usr/share/mime/generic-icons'
IgnorePath '/usr/share/mime/globs' # File extensions
IgnorePath '/usr/share/mime/globs2' # Weighted file extensions?
IgnorePath '/usr/share/mime/icons'
IgnorePath '/usr/share/mime/magic' # Binary magic database
IgnorePath '/usr/share/mime/mime.cache' # Binary
IgnorePath '/usr/share/mime/subclasses'
IgnorePath '/usr/share/mime/treemagic' # Directory magic
IgnorePath '/usr/share/mime/types'
IgnorePath '/usr/share/mime/version'

# fonts
IgnorePath '/etc/fonts/conf.d/*'
IgnorePath '/usr/share/fonts/*'
IgnorePath '/usr/share/glib-2.0/schemas/gschemas.compiled'

# icons
IgnorePath '/usr/share/icons/*/icon-theme.cache'

# docker
IgnorePath '/var/lib/docker*'
IgnorePath '/etc/docker/key.json'

# Artix
IgnorePath '/etc/artix-release'

