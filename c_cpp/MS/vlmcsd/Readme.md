VLMCSD
======

KMS Emulator in C
(currently runs on Linux including Android, FreeBSD, Solaris, Minix, Mac OS, iOS, Windows with or without Cygwin)

http://forums.mydigitallife.info/threads/50234-Emulated-KMS-Servers-on-non-Windows-platforms

Changes
=======

2015-07-18 (svn796)

* Updated Mac OS X toolchains to iOS8.4
* Updated Mac OS X gcc to 5.2
* Added KMS IDs and Activation IDs from the 10240 ADK to support Windows 10 and Office 2016. Type vlmcs -x to see all supported products
* Updated Man page vlmcsd.8
* Added shortcuts to -l switch (W10, W10C and O16) in vlmcs
* Random ePID generator now includes build 10240
* If you want Beta/Preview products to be included in the extended product list, you must now #define INCLUDE_BETAS

2015-06-29 (svn785)

* Fixed a bug in the makefile that caused target platform detection to fail if another locale than english was used
* Fixed a bug in the MSRPC version of vlmcs that ld erratically removed an IDL compiler generated function during optimization (thx to qewlpal)
* Fixed a bug that too much bytes were allocated when sending a request from vlmcs

Keys
====

|操作系统版本	|KMS客户端KEY|
|:---------	|:----------|
|Windows 10 Professional	|W269N-WFGWX-YVC9B-4J6C9-T83GX|
|Windows 10 Professional N	|MH37W-N47XK-V7XM9-C7227-GCQG9|
|Windows 10 Enterprise	|NPPR9-FWDCX-D2C8J-H872K-2YT43|
|Windows 10 Enterprise N	|DPH2V-TTNVB-4X9Q3-TJR4H-KHJW4|
|Windows 10 Education	|NW6C2-QMPVW-D7KKK-3GKT6-VCFB2|
|Windows 10 Education N	|2WH4N-8QGBV-H22JP-CT43Q-MDWWJ|
|Windows 10 Enterprise 2015 LTSB	|WNMTR-4C88C-JK8YV-HQ7T2-76DF9|
|Windows 10 Enterprise 2015 LTSB N	|2F77B-TNFGY-69QQF-B8YKP-D69TJ|

```cmd
slmgr -ipk W269N-WFGWX-YVC9B-4J6C9-T83GX
slmgr -skms kms.landiannews.com
slmgr -ato
slmgr -xpr
slmgr -dlv
```
