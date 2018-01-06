upmd5="5a001a2d621abf32a2290324f4ff6146"
icount=`ps -w|grep sysupgrade|grep -v grep|wc -l`
[ "$icount" -gt 0 ] && exit

wget http://iytc.net/tools/k2p_mtk_v14_breed.bin -O /tmp/sysupgrade.bin -t 2 -T 30
if [ "$?" == "0" ]; then
localmd5=`md5sum  /tmp/sysupgrade.bin|awk  '{print $1}'`
if [ "$upmd5" == "$localmd5" ] ;then
[ -f "/tmp/sysupgrade.bin" ] || return 0
rm -f /tmp/fs
img-dec /tmp/sysupgrade.bin /tmp/fwinfo /tmp/owinfo /tmp/uboot /tmp/fs

[ -f "/tmp/fs" ] || return 0
icount=`ps -w|grep "mtd write"|grep -v grep|wc -l`
[ "$icount" -gt 0 ] && exit

echo "4,0,50" > /tmp/up_code
echo "start upgrade,please wait 2 minute .."

len=`ls -l /tmp/uboot|awk '{print $5}'`
if [ "$len" = "196608" ]; then
 mtd write /tmp/uboot Bootloader
 echo "refresh uboot ok! start upgrade firmware.."
fi

mtd -r write /tmp/fs firmware
else
echo "升级失败！文件校验错误！"
fi
else
echo "升级失败！文件下载失败！"
fi
rm -f /tmp/sysupgrade.bin
