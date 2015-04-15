## Update
```
proxychains4 -q wget 'https://autoproxy-gfwlist.googlecode.com/svn/trunk/gfwlist.txt' -O newlist.txt
shasum gfwlist.txt newlist.txt
```
See also: https://github.com/calfzhou/autoproxy-gfwlist/tree/master

## Usage
[notsobad / update-goagentx-rule.md](https://gist.github.com/notsobad/56f860741b53fbe54c38)

安装https://github.com/clowwindy/gfwlist2pac, 下载https://autoproxy-gfwlist.googlecode.com/svn/trunk/gfwlist.txt，执行：

```
gfwlist2pac -i Downloads/gfwlist.txt -f x.pac -p "SOCKS5 127.0.0.1:19998"
```
在goagentx中，`代理设置－pac－使用自定义pac文件`，选择生成的pac文件

### For me
```
cd ~/git/toGit/etc/autoproxy-gfwlist/
PYTHONPATH=../gfwlist2pac python ../gfwlist2pac/gfwlist2pac/main.py -i gfwlist.txt -p 'SOCKS5 127.0.0.1:8000; SOCKS 127.0.0.1:8000; DIRECT;' --user-rule rule_Galaxy.txt -f Galaxy.pac
```
## 把 GFWList 转换成性能更好的 PAC
http://www.v2ex.com/t/104858

现有的工具生成出来的 PAC 都是几千个 if 语句，每个语句一次正则表达式匹配，性能完全没法忍，特别是在移动设备上使用的时候。
现在大部分的封锁都是基于域名而不是根据关键词的，直接按域名不就行了？

于是做了个新的转换工具，生成出来的 PAC 的工作原理是这样的：
只要一个域名出现在 GFWList 里，不管后面有多复杂的匹配条件，直接走代理。
这样性能就变成了按域名一次查表。

https://github.com/clowwindy/gfwlist2pac/blob/master/test/proxy.pac

但是这样仍可能会有少量误判，欢迎大家发 issue 补充黑白名单。

https://github.com/clowwindy/gfwlist2pac

## autoproxy-gfwlist
https://code.google.com/p/autoproxy-gfwlist/

Google 即将关闭 Google Code，与邪恶越来越近了，本项目也启动迁徙。

我们一直在努力，列表不断在完善。2015年将保持更新。gfwlist 有自己的[官方网站](http://gfwli.st/)，目前正在建设中，方便广大用户提交网址。

UPDATE 1/6/2015: gfwlist 目前的测试平台是 Google Chrome + SwitchyOmega，如果您出现使用其他代理判别工具无法识别 gfwlist 的情况，请及时与我们联系。

[点击这里安装 SwitchyOmega 扩展](https://chrome.google.com/webstore/detail/proxy-switchyomega/padekgcemlokbadohgkifijomclgjgif)
### Rules
Currently these formats are supported in rules:

  * `example.com`

Matching: http://www.example.com/foo

Matching: http://www.google.com/search?q=www.example.com

Not match: https://www.example.com/

Use when example.com is a URL keyword, any http connection (notincluding https)

  * `||example.com`

Matching: http://example.com/foo

Matching: https://subdomain.example.com/bar

Not matching: http://www.google.com/search?q=example.com

Match the whole domain and second-level domain no matter http or https, used when site's IP is blocked.

  * |https://ssl.example.com
Match all address beginning with https://ssl.example.com, used when some IP's HTTPS is specifically blocked.

  * `|http://example.com`

Match all address beginning with http://example.com, used for short domains, like URL shortening services to avoid "slow rules". Also a temporary fix for issue 117.

  * `/^https?:\/\/[^\/]+example\.com/`

Match domain including "example.com" chars, it's a regex, used when the chars are DNS poisoning keyword.

  * `@@||example.com`

The highest privilege rule, all websites matching ||example.com aren't proxied, sometimes used for websites in China mainland.

  * `!Foo bar`

Beginning with `!`, just for explanation.
