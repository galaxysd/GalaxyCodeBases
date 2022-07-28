#!/bin/sh

if [ ! -f /etc/hostname ]; then
	echo 'GalaxyVM' > /etc/hostname
fi

if [ -f /etc/hostname ]; then
	echo -ne '\n# This host address\n127.0.0.1\t' >> /etc/hosts
	cat /etc/hostname >> /etc/hosts
	cp --reflink=auto -Pp /etc/hosts /etc/hosts.init
fi

rm -f $0
