# BASE
a fast and accurate de novo genome assembler for longer NGS reads

`bin/base` is too big (1.6M), thus removed.

from <https://github.com/dhlbh/BASE>, writing by _Binghang Liu_, who once worked in BGI and now works in ICarbonX([深圳碳云智能科技有限公司](http://www.icarbonx.com/)).

Thesis: [Large genome de novo assembly with bi-directional BWT](https://hub.hku.hk/handle/10722/225227) [by](https://hub.hku.hk/bitstream/10722/225227/1/FullText.pdf) Liu, Binghang (劉兵行) @ HKU.

## Comments

BWT 是一种可快速查询子串的字符串索引，从而使 Reads overlap 关系可以 on-the-fly 地构建欧拉图。Seed 得以向两侧实时延伸到分叉为止。

相比 DBG ，which 只能设置固定的 overlap 长度，而非 BWT 在查询中只限定下限。所以对 repeats ，BWT 实现的及时查询具有明显优势。
实际数据表明，当读长大于 100 bp 时，优势可以体现出来。

DBG 的经典实现是 [Velvet](https://github.com/dzerbino/velvet)，[SOAPdenovo](https://github.com/aquaskyline/SOAPdenovo2) 属于在其基础上的简化以压低内存占用。据刘兵行在 2016/5/21 说，DBG 相关概念和实现它那里都有。
