# Surge

Surge is a web developer tool and proxy utility for iOS 9. This app is designed for developers and may require some level of professional knowledge to use.

* <http://surge.run/manual/>
* Trello: <https://trello.com/b/qy4sDvxg/surge>

## SSEncrypt.module of shadowsocks

<http://surge.run/config-example/ss-module.conf>

https://github.com/galaxysd/GalaxyCodeBases/raw/master/etc/Surge/SSEncrypt.module

### Sources

* <https://twitter.com/lovegoodbest/status/649524156589969408>

存档备用。<https://www.dropbox.com/s/v94zjlx2w2yd4za/SSEncrypt.module?dl=0>  
SHA1：400603009de8356d66631b166e5b95b95cdf0cf6

https://github.com/ky0ncheng/gfvvlist/tree/master/surge

* [如何配置surge在iOS中使用shadowsocks科学上网](http://ideafoc.us/2015/10/%E5%A6%82%E4%BD%95%E9%85%8D%E7%BD%AEsurge%E5%9C%A8ios%E4%B8%AD%E4%BD%BF%E7%94%A8shadowsocks%E7%A7%91%E5%AD%A6%E4%B8%8A%E7%BD%91/)

<https://dl.dropboxusercontent.com/u/356582699/SSEncrypt.module>

* [iOS科学上网神器 Surge+SS配置](http://www.jianshu.com/p/41cdb6f71555)

<https://dl.dropboxusercontent.com/u/760466/SSEncrypt.module>

### conf examples

````ini
[Proxy]
Proxy1 = http,$IP,$PORT,$USERNAME,$PASSWORD
Proxy2 = https,$IP,$PORT,$USERNAME,$PASSWORD
Proxy3 = socks,$IP,$PORT
Proxy4 = custom,$IP,$PORT,$METHOD,$PASSWORD,$MODULE_URL
````

````ini
# 主要用于解决 Google.cn 跳转问题，目前该配置项没有 UI。样例配置如下：
[URL Rewrite]
^http://www.google.cn http://www.google.com
````

每行两个参数，以空格分隔，第一个是匹配用的正则表达式，第二个是替换的内容。正则的输入字符串是完整的 URL（包含 http://），输出的字符串也必须以 http:// 开头。该功能只对非 https 请求有效。
