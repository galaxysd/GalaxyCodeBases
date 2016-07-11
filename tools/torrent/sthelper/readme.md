# 1.准备安装
Linux / OSX /windows 用户请先安装 NodeJs 
可以到 https://nodejs.org/en/download/ 下载安装 

## windows用户 
可以下载

64位 https://nodejs.org/dist/v4.4.0/node-v4.4.0-x64.msi
32位 （小白请直接装这个）https://nodejs.org/dist/v4.4.0/node-v4.4.0-x86.msi

安装后记得重启后继续

## Mac OS X 用户 

可以下载 https://nodejs.org/dist/v4.4.0/node-v4.4.0.pkg


## Linux用户
请尽量从源码编译 Debian stable branch / CentOS 等的包管理器从源中安装的Node版本可能过旧
源码例如 https://nodejs.org/dist/v4.4.0/node-v4.4.0.tar.gz

怎么编译我就不赘述了


# 2.安装sthelper

安装完nodeJs后

打开你的终端
windows 按 windows键（键盘上windows图标那个）+ R   输入cmd 回车

mac在启动器找到 终端 这个程序

Linux 呵呵呵我不说了

输入npm install sthelper -g 即可

mac 和Linux 可能要使用 sudo

# 3.使用

安装完成后直接在终端使用sthelper指令即可

使用方法

创建带crc的种子
make [文件目录] [输出种子地址 默认当前目录 out.torrent]
link [数据文件目录] [目标种子文件] [目标链接目录]