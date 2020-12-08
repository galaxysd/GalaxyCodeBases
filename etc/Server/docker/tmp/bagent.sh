#!/bin/bash
SCRIPT_VERSION="1.1.8"

console_domain=""
console_h_param=""
authorization_check_param="1"

download_url_and_param="-s -L 'http://10.225.10.78/agent/installer?k=cc2f168917d40ea73b2d568619f167b81d7896c0&group=1&protocol=0'"
is_ipv6="0"
# remove ' around url
download_url_and_param=`echo ${download_url_and_param} | sed "y/'/ /"`
retry_download_url_and_param=
deploy_server_url=""
logdir=/var/log/titanagent
installlog=$logdir/install.log
proxy_list_conf=/tmp/install_agent/proxy_list.conf

connect_timeout=60
max_time=3600
socks5_option=
if [[ -f /titan/agent/curl ]]; then
    chmod +x /titan/agent/curl
    curl_real=/titan/agent/curl
else
    curl_real=curl
fi

if [[ "$is_ipv6" == "1" ]]; then
curl_real="$curl_real -g -6"
fi

init_h_param()
{
    if [ "$console_domain" != "" ]; then
        console_h_param="-H $console_domain"
    fi
}

init_log()
{
    if [ ! -d $logdir ]; then
        mkdir -p $logdir
    fi
    chmod 700 $logdir
    touch $installlog
    chmod 600 $installlog
}

install_log()
{
    echo "$1" >> $installlog
}

console_log()
{
    local now=`date +"%F %T"`
    install_log "[$$]"" ""$now"" ""- $1"
    if [ -c /dev/tty ]; then
        echo "[$$]"" ""$now"" ""- $1" > /dev/tty
    fi
}

install_log_info()
{
    install_log "[INFO] - $1"
}

install_log_warn()
{
    install_log "[WARN] - $1"
}

install_log_error()
{
    install_log "[ERROR] - $1"
}

console_log_info()
{
    console_log "[INFO] - $1"
}

console_log_warn()
{
    console_log "[WARN] - $1"
}

console_log_error()
{
    console_log "[ERROR] - $1"
}

check_socks5_hostname()
{
    local output=`${curl_real} --socks5-hostname 2>&1 | grep unknown`
    if [[ ! -z $output ]]; then
        socks5_option="--socks5"
    else
        socks5_option="--socks5-hostname"
        retry_download_url_and_param=`echo ${download_url_and_param} | sed "s/--socks5 /--socks5-hostname /"`
    fi
}

fun_after_install()
{
    sed -i -e "s/pass=.*/pass=/g" $installlog
    sed -i -e "s/--proxy-user .* http/--proxy-user ****** http/g" $installlog
    rm -f $proxy_list_conf
}


get_local_ip_list()
{
    local local_ip_list=`ip addr | grep "inet " | grep -v '127.0.0.1' | awk '{print $2}' | sed -e 's/\/.*//'`
    if [ -z $local_ip_list ]; then
        local_ip_list=`ifconfig | grep 'inet ' | grep -v '127.0.0.1' | awk '{print $2}' | sed -e 's/addr://'`
    fi
    echo ${local_ip_list} | sed 's/ /,/g'
}

get_proxy_list()
{
    console_log_info "Download proxy_list from deploy server ${deploy_server_url}"
    local local_ip_list=`get_local_ip_list`

    rm -f $proxy_list_conf
    ${curl_real} -k -v -s --connect-timeout ${connect_timeout} -m ${max_time} -L ${deploy_server_url}/proxy_list?ip_list=${local_ip_list} -o $proxy_list_conf 1>> $installlog 2>> $installlog

    local ret=$?
    if [ ${ret} -ne 0 ]; then
        console_log_error "Error Code 0x00000002: failed to download proxy list, curl exit code is ${ret}"
        return 1
    fi

    if [ ! -s $proxy_list_conf ]; then
        console_log_error "Error Code 0x00000002: download proxy_list error"
        return 1
    fi

    console_log_info "Download proxy_list successfully"
    return 0
}

do_download_install_agent()
{
    local url_and_param=$1
    ${curl_real} -k -v -s --connect-timeout ${connect_timeout} -m ${max_time} -L ${url_and_param} ${console_h_param} -o /tmp/install_agent/install_agent.sh 1>> $installlog 2>> $installlog

    local ret=$?
    if [ ${ret} -ne 0 ]; then
        console_log_error "Error Code 0x00000003: failed to download intall_agent.sh, curl exit code is ${ret}"
        return 1
    fi

    if [ ! -f /tmp/install_agent/install_agent.sh ]; then
        console_log_error "Error Code 0x00000003: download intall_agent.sh error"
        return 1
    fi

    console_log_info "Download install_agent.sh successfully"
    return 0
}

download_install_agent()
{
    if [ -f $proxy_list_conf ]; then
        while read line
        do
            local proxy_ip=`echo ${line} | cut -d: -f1`
            local proxy_port=`echo ${line} | cut -d: -f2`
            local proxy_user=`echo ${line} | cut -d: -f3`
            local proxy_passwd=`echo ${line} | cut -d: -f4`

            local proxy_info
            local log_proxy_info
            if [ ${proxy_ip} != "direct" ]; then
                proxy_info="${socks5_option} ${proxy_ip}:${proxy_port} --proxy-user ${proxy_user}:${proxy_passwd}"
                log_proxy_info="${socks5_option} ${proxy_ip}:${proxy_port} --proxy-user ${proxy_user}:******"
            else
                proxy_info=
            fi

            install_log "download install_agent.sh, param ${download_url_and_param} ${log_proxy_info}"
            if do_download_install_agent "${download_url_and_param} ${proxy_info}"; then
                if ! sed -i -e "s/^proxy_ip=.*/proxy_ip=${proxy_ip}/" -e "s/^proxy_port=.*/proxy_port=${proxy_port}/" -e "s/^proxy_user=.*/proxy_user=${proxy_user}/" -e "s:^proxy_pwd=.*:proxy_pwd=${proxy_passwd}:" -e "s#^deploy_server_url=.*#deploy_server_url=${deploy_server_url}#" /tmp/install_agent/install_agent.sh; then
                    install_log_error "failed to write proxy info and deploy server url to install_agent.sh"
                    return 1
                fi
                return 0
            fi
        done < $proxy_list_conf
    else
        install_log "download install_agent.sh, param ${download_url_and_param}"
        do_download_install_agent "${download_url_and_param}"
        if [[ $? -ne 0 ]] && [[ ! -z $retry_download_url_and_param ]]; then
            install_log_warn "failed to download install_agent.sh with --socks5 option, try to use --socks5-hostname option"
            install_log "download install_agent.sh, param ${retry_download_url_and_param}"
            do_download_install_agent "${retry_download_url_and_param}"
        fi
        return $?
    fi

    return 1
}

install_agent()
{
    cd /tmp/install_agent/
    sed -i '/^domain_map_table/s/\"/\\\"/g' /tmp/install_agent/install_agent.sh
    bash ./install_agent.sh install 1>> $installlog 2>> $installlog
}

check_environment()
{
    if [ `id -u` -ne 0 ]; then
        console_log_error "Titan Agent installation need running as root."
        exit 1
    fi

    uname -m | grep 64 > /dev/null
    if [ $? -ne 0 ]; then
        console_log_info "Titan Agent Running on a 32-bit platform."
        # do something...
    else
        console_log_info "Titan Agent Running on a 64-bit platform."
        # do something...
    fi

    if [ "$authorization_check_param" -ne "1" ]; then
        console_log_error "Error Code 0x00000001: license not enough, Can not install"
        exit 1
    fi
}

init_log
check_environment

rm -rf /tmp/install_agent
mkdir /tmp/install_agent
chmod 700 /tmp/install_agent

init_h_param
check_socks5_hostname
console_log_info "Titan Agent $SCRIPT_VERSION installation start"
if [ -n "${deploy_server_url}" ]; then
    if ! get_proxy_list; then
        console_log_error "Error: An error occurred when installing Titan Agent. Please go to /var/log/titanagent/install.log for more information."
        exit 1
    fi
fi

if download_install_agent; then
    install_agent
else
    console_log_error "Error: An error occurred when installing Titan Agent. Please go to /var/log/titanagent/install.log for more information."
fi

fun_after_install
