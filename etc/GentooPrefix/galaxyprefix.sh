#!/bin/env bash

# tar -Pczvf galaxyprefix.tgz /etc/make.conf /etc/portage /opt/gentoo/etc/make.conf /opt/gentoo/etc/portage /opt/gentoo/etc/eixrc
#lr--r--r-- 1 root +Administrators 25 Jan 21 04:32 /etc/make.conf -> /opt/gentoo/etc/make.conf
#lr--r--r-- 1 root +Administrators 23 Jan 21 04:33 /etc/portage -> /opt/gentoo/etc/portage

emerge -a dev-python/pysqlite
mv /opt/gentoo/var/cache/edb/dep/opt /opt/gentoo/var/cache/edb/dep/opt0
emerge --meta

tar -Pxzvf galaxyprefix.tgz

