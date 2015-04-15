## Update
```
proxychains4 -q wget 'https://autoproxy-gfwlist.googlecode.com/svn/trunk/gfwlist.txt' -O newlist.txt
shasum gfwlist.txt newlist.txt
```
See also: https://github.com/calfzhou/autoproxy-gfwlist/tree/master

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

### Rules
Currently these formats are supported in rules:

  * `example.com`

Matching: http://www.example.com/foo

Matching: http://www.google.com/search?q=www.example.com

Not match: https://www.example.com/

Use when example.com is a URL keyword, any http connection (notincluding https)

  * ||example.com

Matching: http://example.com/foo
Matching: https://subdomain.example.com/bar
Not matching: http://www.google.com/search?q=example.com
Match the whole domain and second-level domain no matter http or https, used when site's IP is blocked.

  * |https://ssl.example.com
Match all address beginning with https://ssl.example.com, used when some IP's HTTPS is specifically blocked.

  * |http://example.com
Match all address beginning with http://example.com, used for short domains, like URL shortening services to avoid "slow rules". Also a temporary fix for issue 117.

  * /^https?:\/\/[^\/]+example\.com/
Match domain including "example.com" chars, it's a regex, used when the chars are DNS poisoning keyword.

  * @@||example.com
The highest privilege rule, all websites matching ||example.com aren't proxied, sometimes used for websites in China mainland.

  * !Foo bar
Beginning with !, just for explanation.
