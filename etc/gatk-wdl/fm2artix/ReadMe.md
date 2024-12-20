## Server

```bash
cd /share/dat0/cromwell
#cp -a ../test/bin .
#ln -s ../fq .
cromwell server
```

## Client

```bash
zip -9r fm2wdl.zip fm2.wdl  tasks/
#cromwell submit --workflow-root /share/dat0/test/ -i fmtest.json /share/dat0/test/fm2.wdl -o cromwell_options_no_user.json -p fm2wdl.zip
cromwell submit -i fmtest.json fm2.wdl -o cromwell_options_no_user.json -p fm2wdl.zip
```

## Status

<http://192.168.2.151:8000/api/workflows/v1/ca03f664-1d02-414d-a566-d496bcb51cf6/timing>

## Install

```bash
pacman -S python-yaml
expac --timefmt='%F %T' '%l %n' | sort -n > expac.lst
yay -S jre11-openjdk imagej
pacman -S adobe-source-han-sans-otc-fonts adobe-source-han-serif-otc-fonts
```

```
2021-07-06 02:21:04 filesystem
2021-07-06 02:21:04 glibc
2021-07-06 02:21:04 iana-etc
2021-07-06 02:21:04 linux-api-headers
2021-07-06 02:21:04 tzdata
2021-07-06 02:21:05 acl
2021-07-06 02:21:05 attr
2021-07-06 02:21:05 audit
2021-07-06 02:21:05 bash
2021-07-06 02:21:05 e2fsprogs
2021-07-06 02:21:05 gcc-libs
2021-07-06 02:21:05 gmp
2021-07-06 02:21:05 keyutils
2021-07-06 02:21:05 krb5
2021-07-06 02:21:05 libcap-ng
2021-07-06 02:21:05 libldap
2021-07-06 02:21:05 libsasl
2021-07-06 02:21:05 libtirpc
2021-07-06 02:21:05 libxcrypt
2021-07-06 02:21:05 ncurses
2021-07-06 02:21:05 openssl
2021-07-06 02:21:05 pam
2021-07-06 02:21:05 pambase
2021-07-06 02:21:05 readline
2021-07-06 02:21:05 util-linux-libs
2021-07-06 02:21:06 bzip2
2021-07-06 02:21:06 coreutils
2021-07-06 02:21:06 file
2021-07-06 02:21:06 findutils
2021-07-06 02:21:06 gawk
2021-07-06 02:21:06 gettext
2021-07-06 02:21:06 glib2
2021-07-06 02:21:06 grep
2021-07-06 02:21:06 hwids
2021-07-06 02:21:06 icu
2021-07-06 02:21:06 kmod
2021-07-06 02:21:06 libcroco
2021-07-06 02:21:06 libelogind
2021-07-06 02:21:06 libeudev
2021-07-06 02:21:06 libffi
2021-07-06 02:21:06 libseccomp
2021-07-06 02:21:06 libunistring
2021-07-06 02:21:06 libxml2
2021-07-06 02:21:06 lz4
2021-07-06 02:21:06 mpfr
2021-07-06 02:21:06 pciutils
2021-07-06 02:21:06 pcre
2021-07-06 02:21:06 procps-ng
2021-07-06 02:21:06 psmisc
2021-07-06 02:21:06 sed
2021-07-06 02:21:06 tar
2021-07-06 02:21:06 xz
2021-07-06 02:21:06 zlib
2021-07-06 02:21:06 zstd
2021-07-06 02:21:07 ca-certificates
2021-07-06 02:21:07 ca-certificates-mozilla
2021-07-06 02:21:07 ca-certificates-utils
2021-07-06 02:21:07 curl
2021-07-06 02:21:07 expat
2021-07-06 02:21:07 gzip
2021-07-06 02:21:07 less
2021-07-06 02:21:07 libarchive
2021-07-06 02:21:07 libassuan
2021-07-06 02:21:07 libgcrypt
2021-07-06 02:21:07 libgpg-error
2021-07-06 02:21:07 libidn2
2021-07-06 02:21:07 libksba
2021-07-06 02:21:07 libnghttp2
2021-07-06 02:21:07 libp11-kit
2021-07-06 02:21:07 libpsl
2021-07-06 02:21:07 libsecret
2021-07-06 02:21:07 libssh2
2021-07-06 02:21:07 libtasn1
2021-07-06 02:21:07 licenses
2021-07-06 02:21:07 nettle
2021-07-06 02:21:07 npth
2021-07-06 02:21:07 p11-kit
2021-07-06 02:21:07 pcre2
2021-07-06 02:21:07 pinentry
2021-07-06 02:21:07 shadow
2021-07-06 02:21:07 util-linux
2021-07-06 02:21:08 artix-keyring
2021-07-06 02:21:08 dbus
2021-07-06 02:21:08 eudev
2021-07-06 02:21:08 gnupg
2021-07-06 02:21:08 gnutls
2021-07-06 02:21:08 gpgme
2021-07-06 02:21:08 kbd
2021-07-06 02:21:08 kexec-tools
2021-07-06 02:21:08 libx11
2021-07-06 02:21:08 libxcb
2021-07-06 02:21:08 pacman
2021-07-06 02:21:08 sqlite
2021-07-06 02:21:08 xcb-proto
2021-07-06 02:21:08 xorgproto
2021-07-06 02:21:09 base
2021-07-06 02:21:09 db
2021-07-06 02:21:09 dbus-openrc
2021-07-06 02:21:09 diffutils
2021-07-06 02:21:09 elogind
2021-07-06 02:21:09 elogind-openrc
2021-07-06 02:21:09 gdbm
2021-07-06 02:21:09 inetutils
2021-07-06 02:21:09 iproute2
2021-07-06 02:21:09 iptables
2021-07-06 02:21:09 iputils
2021-07-06 02:21:09 libelf
2021-07-06 02:21:09 libmnl
2021-07-06 02:21:09 libnetfilter_conntrack
2021-07-06 02:21:09 libnfnetlink
2021-07-06 02:21:09 libnftnl
2021-07-06 02:21:09 libnl
2021-07-06 02:21:09 libpcap
2021-07-06 02:21:09 m4
2021-07-06 02:21:09 netifrc
2021-07-06 02:21:09 openrc
2021-07-06 02:21:10 autoconf
2021-07-06 02:21:10 automake
2021-07-06 02:21:10 binutils
2021-07-06 02:21:10 bison
2021-07-06 02:21:10 elfutils
2021-07-06 02:21:10 fakeroot
2021-07-06 02:21:10 flex
2021-07-06 02:21:10 libmpc
2021-07-06 02:21:10 perl
2021-07-06 02:21:11 gcc
2021-07-06 02:21:11 groff
2021-07-06 02:21:11 libtool
2021-07-06 02:21:12 gc
2021-07-06 02:21:12 guile
2021-07-06 02:21:12 make
2021-07-06 02:21:12 patch
2021-07-06 02:21:12 pkgconf
2021-07-06 02:21:12 sudo
2021-07-06 02:21:12 which
2021-07-06 14:03:08 gptfdisk
2021-07-06 18:06:09 mkinitcpio-busybox
2021-07-06 18:06:11 linux-firmware
2021-07-06 18:14:03 device-mapper
2021-07-06 18:14:03 dnssec-anchors
2021-07-06 18:14:03 gpm
2021-07-06 18:14:03 grub
2021-07-06 18:14:03 ldns
2021-07-06 18:14:03 libedit
2021-07-06 18:14:03 openssh
2021-07-06 18:14:03 openssh-openrc
2021-07-06 18:14:03 vim
2021-07-06 18:14:03 vim-runtime
2021-07-06 18:14:03 wpa_supplicant
2021-07-06 18:14:03 wpa_supplicant-openrc
2021-07-06 18:18:49 libnsl
2021-07-06 18:18:50 dosfstools
2021-07-06 18:18:50 efibootmgr
2021-07-06 18:18:50 efivar
2021-07-06 18:18:50 freetype2
2021-07-06 18:18:50 fuse-common
2021-07-06 18:18:50 fuse2
2021-07-06 18:18:50 graphite
2021-07-06 18:18:50 libburn
2021-07-06 18:18:50 libisoburn
2021-07-06 18:18:50 libisofs
2021-07-06 18:18:50 libpng
2021-07-06 18:18:50 os-prober
2021-07-06 18:18:50 popt
2021-07-06 18:24:01 pahole
2021-07-06 18:24:03 btrfs-progs
2021-07-06 18:24:03 dkms
2021-07-06 18:24:03 linux-headers
2021-07-06 18:24:03 lzo
2021-07-06 18:24:43 intel-ucode
2021-07-06 18:24:43 linux
2021-07-06 18:24:43 memtest86+
2021-07-06 18:27:09 acpid
2021-07-06 18:27:09 acpid-openrc
2021-07-06 18:27:09 cronie
2021-07-06 18:27:09 cronie-openrc
2021-07-06 18:27:09 dhcpcd-openrc
2021-07-06 18:27:09 haveged
2021-07-06 18:27:09 haveged-openrc
2021-07-06 18:27:09 json-c
2021-07-06 18:27:09 libnet
2021-07-06 18:27:09 logrotate
2021-07-06 18:27:09 ntp
2021-07-06 18:27:09 ntp-openrc
2021-07-06 18:27:09 openresolv
2021-07-06 18:27:09 run-parts
2021-07-06 18:27:09 syslog-ng-openrc
2021-07-07 15:39:41 python
2021-07-07 15:39:41 texinfo
2021-07-09 14:33:50 dhcpcd
2021-07-09 14:33:50 esysusers
2021-07-09 14:33:50 etmpfiles
2021-07-09 14:33:50 libcap
2021-07-09 14:33:50 mkinitcpio
2021-07-12 10:22:04 artix-mirrorlist
2021-07-12 10:22:05 harfbuzz
2021-07-12 10:22:05 syslog-ng
2021-07-12 11:16:35 postgresql-libs
2021-07-12 11:16:36 libxslt
2021-07-12 11:16:36 llvm-libs
2021-07-12 11:16:36 postgresql
2021-07-12 11:16:36 postgresql-openrc
2021-07-12 14:21:19 expect
2021-07-12 14:21:19 nspr
2021-07-12 14:21:19 nss
2021-07-12 14:21:19 oath-toolkit
2021-07-12 14:21:19 tcl
2021-07-12 14:21:19 xmlsec
2021-07-12 15:10:43 libpipeline
2021-07-12 15:10:43 man-db
2021-07-12 15:11:04 rsync
2021-07-12 15:11:04 xxhash
2021-07-12 15:42:14 git
2021-07-12 15:42:14 perl-error
2021-07-12 15:42:14 perl-mailtools
2021-07-12 15:42:14 perl-timedate
2021-07-12 15:52:09 go
2021-07-12 15:56:48 yay
2021-07-12 15:58:41 blas
2021-07-12 15:58:41 cblas
2021-07-12 15:58:41 gsl
2021-07-12 15:58:41 lapack
2021-07-12 15:58:41 python-cycler
2021-07-12 15:58:41 python-dateutil
2021-07-12 15:58:41 python-kiwisolver
2021-07-12 15:58:41 python-six
2021-07-12 15:58:42 fribidi
2021-07-12 15:58:42 lcms2
2021-07-12 15:58:42 libimagequant
2021-07-12 15:58:42 libjpeg-turbo
2021-07-12 15:58:42 libraqm
2021-07-12 15:58:42 libtiff
2021-07-12 15:58:42 openjpeg2
2021-07-12 15:58:42 python-matplotlib
2021-07-12 15:58:42 python-numpy
2021-07-12 15:58:42 python-pyparsing
2021-07-12 15:58:42 qhull
2021-07-12 16:11:40 htslib
2021-07-12 16:12:49 samtools
2021-07-12 16:13:40 bcftools
2021-07-12 16:17:00 java-environment-common
2021-07-12 16:17:00 java-runtime-common
2021-07-12 16:17:00 jdk8-openjdk
2021-07-12 16:17:00 jre8-openjdk
2021-07-12 16:17:00 jre8-openjdk-headless
2021-07-12 16:17:00 libxt
2021-07-12 16:17:00 xorg-xprop
2021-07-12 16:34:42 gatk-bin
2021-07-12 16:39:27 cairo
2021-07-12 16:39:27 fontconfig
2021-07-12 16:39:27 libdatrie
2021-07-12 16:39:27 libthai
2021-07-12 16:39:27 pango
2021-07-12 16:39:27 r
2021-07-12 16:39:27 unzip
2021-07-12 16:39:27 zip
2021-07-12 16:39:42 littler
2021-07-12 16:42:52 imagej
2021-07-12 16:46:15 bwa-git
2021-07-12 16:48:50 wget
2021-07-12 16:54:28 gradle
2021-07-12 17:03:27 maven
2021-07-12 17:25:57 r-digest
2021-07-12 17:26:04 r-glue
2021-07-12 17:26:09 r-gtable
2021-07-12 17:26:33 r-isoband
2021-07-12 17:26:45 r-rlang
2021-07-12 17:27:09 r-farver
2021-07-12 17:27:14 r-labeling
2021-07-12 17:27:18 r-lifecycle
2021-07-12 17:27:29 r-colorspace
2021-07-12 17:27:33 r-munsell
2021-07-12 17:27:38 r-r6
2021-07-12 17:27:42 r-rcolorbrewer
2021-07-12 17:27:46 r-viridislite
2021-07-12 17:27:53 r-scales
2021-07-12 17:27:59 r-ellipsis
2021-07-12 17:28:08 r-fansi
2021-07-12 17:28:14 r-magrittr
2021-07-12 17:28:24 r-cli
2021-07-12 17:28:29 r-crayon
2021-07-12 17:28:38 r-utf8
2021-07-12 17:28:56 r-vctrs
2021-07-12 17:29:02 r-pillar
2021-07-12 17:29:06 r-pkgconfig
2021-07-12 17:29:14 r-tibble
2021-07-12 17:29:19 r-withr
2021-07-12 17:29:36 r-ggplot2
2021-07-12 18:22:27 python2
2021-07-12 18:22:38 cfv
2021-07-13 09:56:00 screen
2021-07-13 10:22:30 hicolor-icon-theme
2021-07-13 10:22:30 libice
2021-07-13 10:22:30 libsm
2021-07-13 10:22:30 libxau
2021-07-13 10:22:30 libxdmcp
2021-07-13 10:22:30 libxext
2021-07-13 10:22:30 libxft
2021-07-13 10:22:30 libxmu
2021-07-13 10:22:30 libxrender
2021-07-13 10:22:30 mtools
2021-07-13 10:22:30 pixman
2021-07-13 10:22:30 xdg-utils
2021-07-13 10:22:30 xorg-xset
2021-07-13 14:10:52 python-appdirs
2021-07-13 14:10:52 python-cachecontrol
2021-07-13 14:10:52 python-cffi
2021-07-13 14:10:52 python-chardet
2021-07-13 14:10:52 python-colorama
2021-07-13 14:10:52 python-contextlib2
2021-07-13 14:10:52 python-cryptography
2021-07-13 14:10:52 python-distlib
2021-07-13 14:10:52 python-distro
2021-07-13 14:10:52 python-html5lib
2021-07-13 14:10:52 python-idna
2021-07-13 14:10:52 python-more-itertools
2021-07-13 14:10:52 python-msgpack
2021-07-13 14:10:52 python-ordered-set
2021-07-13 14:10:52 python-packaging
2021-07-13 14:10:52 python-pep517
2021-07-13 14:10:52 python-pip
2021-07-13 14:10:52 python-ply
2021-07-13 14:10:52 python-progress
2021-07-13 14:10:52 python-pycparser
2021-07-13 14:10:52 python-pyopenssl
2021-07-13 14:10:52 python-requests
2021-07-13 14:10:52 python-resolvelib
2021-07-13 14:10:52 python-retrying
2021-07-13 14:10:52 python-setuptools
2021-07-13 14:10:52 python-toml
2021-07-13 14:10:52 python-urllib3
2021-07-13 14:10:52 python-webencodings
2021-07-13 14:18:23 htop
2021-07-13 14:18:37 libunwind
2021-07-13 14:18:37 lsof
2021-07-13 14:18:37 strace
2021-07-13 14:18:43 lm_sensors
2021-07-13 14:30:52 giflib
2021-07-13 14:30:52 jdk11-openjdk
2021-07-13 14:30:52 jre11-openjdk
2021-07-13 14:30:52 jre11-openjdk-headless
2021-07-13 14:40:00 ttf-dejavu
2021-07-13 14:40:46 fastqc
2021-07-13 14:49:28 python-dataclasses
2021-07-13 14:49:44 cython
2021-07-13 14:51:52 sambamba
2021-07-13 17:17:10 python-pillow
2021-07-13 18:12:05 libyaml
2021-07-13 18:12:05 python-yaml
2021-07-13 18:16:09 expac
```

```bash
$ ls -1 ~/.local/lib/python3.9/site-packages/ |grep info
PyYAML-5.4.1.dist-info
biowdl_input_converter-0.3.0.dist-info
cutadapt-3.4.dist-info
dnaio-0.5.1.dist-info
isal-0.10.0.dist-info
xopen-1.1.0.dist-info
```
