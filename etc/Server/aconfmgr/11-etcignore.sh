# self rm files from 'artix-branding-base'
IgnorePath '/etc/local.d/0PS1.start'
IgnorePath '/etc/local.d/consolefont.start'
IgnorePath '/etc/local.d/mkinitcpio.start'

IgnorePath '/etc/fstab'
IgnorePath '/etc/group*'
IgnorePath '/etc/gshadow*'
IgnorePath '/etc/passwd*'
IgnorePath '/etc/shadow*'
IgnorePath '/etc/xml/catalog'
IgnorePath '/etc/resolv.conf'
IgnorePath '/etc/*.bak'
IgnorePath '/etc/shells'
IgnorePath '__pycache__'
IgnorePath '/etc/.pwd.lock'
IgnorePath '/etc/dconf/db/*'
IgnorePath '/usr/lib/locale/locale-archive'
IgnorePath '/var/lib/sddm/.cache/*'

IgnorePath '/etc/issue'
IgnorePath '/etc/lsb-release'
IgnorePath '/etc/machine-id'
IgnorePath '/etc/NetworkManager/system-connections/*'
IgnorePath '/etc/ssh/ssh_host_*'
IgnorePath '/var/lib/NetworkManager/*'
IgnorePath '/var/lib/random-seed'
IgnorePath '/var/lib/sddm/*'
IgnorePath '/var/lib/upower/*'
IgnorePath '/var/spool/*'
IgnorePath '/var/lib/logrotate.status'
IgnorePath '/etc/cups/subscriptions.conf*'
IgnorePath '/etc/printcap'
IgnorePath '/usr/lib/os-release'
IgnorePath '/etc/machine-id'
IgnorePath '/var/lib/dbus/machine-id'
IgnorePath '/var/lib/syslog-ng/syslog-ng.persist'

# Galaxy
IgnorePath '/.snap*'
IgnorePath '/etc/skel/.bashrc'
