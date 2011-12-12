#!/bin/sh

autossh -M20000 -CNfR 2260:localhost:5321 tun@galaxy.3322.org -p 1233 -i /etc/autossh/tun
autossh -M20002 -CNfR 2261:192.168.8.184:22 tun@galaxy.3322.org -p 1233 -i /etc/autossh/tun
autossh -M20004 -CNfR 2262:192.168.8.185:22 tun@galaxy.3322.org -p 1233 -i /etc/autossh/tun

