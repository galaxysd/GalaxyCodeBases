#!/bin/bash
#
# Author:  Wind4 (puxiaping@gmail.com)
#
# Installs a vlmcsd server for CentOS 6/7
 
# declare
WORK_DIR=/tmp/vlmcsd
 
VL_VERSION="779"
 
# user settings
read -p "Please input vlmcsd SVN version (${VL_VERSION}): " input
VL_VERSION=${input:-$VL_VERSION}
 
# switch work dir
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

# download
yum install wget -y
wget http://files.crsoo.com/packages/vlmcsd/svn${VL_VERSION}/binaries/Linux/intel/static/vlmcsd-x64-musl-static -O ./vlmcsd.run || { echo "Download failed"; exit 1; }
wget https://gist.github.com/Wind4/6ad83d2bea0d3e49a5e3/raw/c011bb4c5992602bbbfe743e57ded0d6869beb88/vlmcsd -O ./vlmcsd.init || { echo "Download failed"; exit 1; }

# config
cp ./vlmcsd.run //usr/bin/vlmcsd
chmod +x /usr/bin/vlmcsd
cp ./vlmcsd.init /etc/init.d/vlmcsd
chmod +x /etc/init.d/vlmcsd
chkconfig --add vlmcsd
service vlmcsd start
 
# clean up
rm -rf ${WORKDIR}