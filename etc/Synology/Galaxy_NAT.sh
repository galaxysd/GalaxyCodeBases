#!/bin/bash

#
# Change this variable to match your private network.
#
PRIVATE_NETWORK="172.99.0.0/16"

#
# Change this variable to match your public interface - either eth0 or eth1
#
PUBLIC_INTERFACE="eth0"

#
# Set PATH to find iptables
#
#PATH=/sbin:/bin:/usr/sbin:/usr/bin:/usr/syno/sbin:/usr/syno/bin

#
# Module list where KERNEL_MODULES_NAT are defined.
#
#IPTABLES_MODULE_LIST="/usr/syno/etc/iptables_modules_list"
IPTABLES_MODULE_LIST="/usr/syno/etc.defaults/iptables_modules_list"
#source "${IPTABLES_MODULE_LIST}"

#
# My service name - let's make sure we don't conflict with synology
#
SERVICE="Galaxy_NAT"

#
# iptable binary
#
IPTABLES="/sbin/iptables"

BIN_IPTABLESTOOL="/usr/syno/bin/iptablestool"
BIN_SYNOMODULETOOL="/usr/syno/bin/synomoduletool"

reverse_modules() {
	local modules=$1
	local mod
	local ret=""
	for mod in $modules; do
	    ret="$mod $ret"
	done
	echo $ret
}

# Based on /volume1/@appstore/VPNCenter/scripts/accel-pppd.sh
NAT_Mod=""
if [ -f "${IPTABLES_MODULE_LIST}" ]; then
	source ${IPTABLES_MODULE_LIST}

	for mod in $KERNEL_MODULES_CORE; do
		if [ -e "/lib/modules/$mod" ]; then
			NAT_Mod="${NAT_Mod} ${mod}"
		fi
	done
	for mod in $KERNEL_MODULES_COMMON; do
		if [ -e "/lib/modules/$mod" ]; then
			NAT_Mod="${NAT_Mod} ${mod}"
		fi
	done
	for mod in $KERNEL_MODULES_NAT; do
		if [ -e "/lib/modules/$mod" ]; then
			NAT_Mod="${NAT_Mod} ${mod}"
		fi
	done
else
	echo >&2 "[x]Cannot find ${IPTABLES_MODULE_LIST} !"
	exit 1
fi

start() {
    #
    # Log execution time
    #
    date
    #
    # Make sure packet forwarding is enabled.
    # 'sysctl -w net.ipv4.ip_forward=1' does not work for me
    #
    echo 1 > /proc/sys/net/ipv4/ip_forward
    #
    # Count the number of modules so that we can verify if the module
    # insertion was successful. We replace whitespaces with newlines
    # and count lines.
    #
    MODULE_COUNT=( ${NAT_Mod} )
	MODULE_COUNT=${#MODULE_COUNT[@]}
    #
    # Load the kernel modules necessary for NAT
    #
	echo -n "Starting ${SERVICE}: "
	if [ -x ${BIN_SYNOMODULETOOL} ]; then
		$BIN_SYNOMODULETOOL --insmod $SERVICE ${NAT_Mod}
		RV=$?
	elif [ -x ${BIN_IPTABLESTOOL} ]; then
		$BIN_IPTABLESTOOL --insmod $SERVICE ${NAT_Mod}
		RV=$?
	fi
    #
    # $BIN_SYNOMODULETOOL returns the number of loaded modules as return value
    #
    [[ "${RV}" == "${MODULE_COUNT}" ]] || {
            echo >&2 "Error: Modules were not loaded (${RV},${MODULE_COUNT}). The following command failed:"
            echo >&2 "${BIN_SYNOMODULETOOL}" --insmod "${SERVICE}" ${NAT_Mod}
            exit 1
    }
    #
    # Turn on NAT.
    #
	${IPTABLES} -t nat -F
    "${IPTABLES}" -t nat -A POSTROUTING -s "${PRIVATE_NETWORK}" -j MASQUERADE -o "${PUBLIC_INTERFACE}"
    RV=$?
    [[ "${RV}" == "0" ]] || {
            echo >&2 "Error: MASQUERADE rules could not be added. The following command failed:"
            echo >&2 "${IPTABLES}" -t nat -A POSTROUTING -s "${PRIVATE_NETWORK}" -j MASQUERADE -o "${PUBLIC_INTERFACE}"
	        exit 1
	}
	#
	# Port Forwarding 
	${IPTABLES} -t nat -A PREROUTING -p tcp -i eth0 --dport 222 -j DNAT --to-destination 172.99.3.3:22
    #
    # Log current nat table
    #
    ${IPTABLES} -L -v -t nat --line-numbers
	echo "done."
}
stop() {
	local modules=`reverse_modules "${NAT_Mod}"`
	echo -n "Shutting down ${SERVICE}: "
	# https://www.digitalocean.com/community/tutorials/how-to-list-and-delete-iptables-firewall-rules
	${IPTABLES} -P INPUT ACCEPT
	${IPTABLES} -P FORWARD ACCEPT
	${IPTABLES} -P OUTPUT ACCEPT
	${IPTABLES} -t nat -F
	${IPTABLES} -t mangle -F
	${IPTABLES} -F
	${IPTABLES} -X
	#echo 0 > /proc/sys/net/ipv4/ip_forward
	if [ -x ${BIN_SYNOMODULETOOL} ]; then
		$BIN_SYNOMODULETOOL --rmmod $SERVICE $modules
	elif [ -x ${BIN_IPTABLESTOOL} ]; then
		$BIN_IPTABLESTOOL --rmmod $SERVICE $modules
	fi
	echo "done."
}

case "$1" in
	start)
		start
	;;
	stop)
		stop
	;;
	restart|reload)
		stop
		start
	;;
	*)
		echo "Usage: $0 {start|stop|restart}"
		exit 1
	;;
esac
exit $?
