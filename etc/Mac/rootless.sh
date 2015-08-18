#!/bin/bash
ME=`basename $0`
IS=$(nvram boot-args | awk '{ print $2 }' | grep rootless | cut -d= -f2)
if [ $# -lt 1 ]; then
  if [ "$IS" == "1" ]; then
    echo "Rootless mode currently active. Use $ME 0 to deactivate."
  else
    echo "rootless mode currently inactive. Use $ME 1 to activate."
  fi
  exit
fi
if [ "$1" == "0" ]; then
  if [ "$IS" == "1" ]; then
    sudo nvram boot-args="rootless=0";
    echo "Rootless mode deactivated. Reboot to apply.";
  else
    echo "Rootless mode already inactive. Nothing changed.";
  fi
else
  if [ "$IS" == "0" ]; then
    sudo nvram boot-args="rootless=1";
    echo "Rootless mode activated. Reboot to apply.";
  else
    echo "Rootless mode already active. Nothing changed.";
  fi
fi
