## 如何在WORD中设置标题1与标题2编号样式不一样
http://blog.csdn.net/kevinhg/article/details/6313042

论文要求每一章的大标题编号要求使用中文数字，而其它更小的标题用阿拉伯数字编号，即

```
第一章 标题1 
1.1 标题2 
1.1.1 标题3
```

找上届的师兄借了个模板过来，发现里面的大标题是手动添加的，小标题是全是通过域代码添加的。

当然域代码或者工都可以实现，但维护起来太不方便，我相信WORD这么一个成熟的产品，肯定有办法完成我们想要的功能，因为这功能要求并不过份，于时到网上查资料，GOOGLE半天也没出来有用的结果。
不过有一个帖子提到使用自动编号的高级选项里面的“正规形式编号”的选项，但也没有给出正确的方法，我试了一下，在大标题样式标识为中文数字时，一勾上“正规形式编号”，发现编号格式里自动将中文数字改为了阿拉伯数字，而且编号形式列表被禁用。
于是我想，如果小标题采用“正规形式编号”，中文数字不也成了阿拉伯数字了，不就达到目的了？果然，在一级标题设为中文数字的情况下，把二级标题的“正规形式编号”勾上，发现果然中文数字了，OK！

总结：

1. 设置一级标题编号为中文数字格式。 
2. 设置二级及以下标题为多级编号，在标题上点右键，选择“项目符号和编号”，选择“多级符号”，点击“自定义”，点击“高级”，勾上“正规形式编号”。

## Word下实现“图一.1”变“图1.1”

做了“一”到“十”的嵌套`if`:

```
{IF {STYLEREF 1 \s} = "一" 1 {IF {STYLEREF 1 \s} = "二" 2 {IF {STYLEREF 1 \s} = "三" 3 {IF {STYLEREF 1 \s} = "四" 4 {IF {STYLEREF 1 \s} = "五" 5 {IF {STYLEREF 1 \s} = "六" 6 {IF {STYLEREF 1 \s} = "七" 7 {IF {STYLEREF 1 \s} = "八" 8 {IF {STYLEREF 1 \s} = "九" 9 {IF {STYLEREF 1 \s} = "十" 10 0 }}}}}}}}} \* MERGEFORMAT}
```
结尾的`\* MERGEFORMAT`应该去掉，否则设格式有问题，大概。

2008年，在网上出现了其他方法，[是用汉字日期的](http://club.excelhome.net/thread-312335-1-1.html)，但2011中文Mac版不支持汉字日期转换。
帖子如下：

#### [原创.春日偶成]域代码将题注“图一-1”变成“图1-1”

一、源起：

1. `VBA.IsDate()`，在`IsDate()`函数中，我们可以发现Word似乎可以判断诸如“一九八一年十月一日”这样的日期数据（`VBA.IsDate("一九八一年十月一日")=True`.
2. 在ASK域中，Word可以根据输入的日期值，转换为日期数据，由此联想到了SET域，结果发现域代码:`{ SET myBK "一九一一年一月一日" }{ myBK \@ "D" }`值为“1”，当然，其间我测试了N遍，从年到月到日，最终确定使用日的范围更广一些（极限值为31，通常对于Word写作而言，三十一章（标题）基本适用了）。

二、StringNumber

也许大家知道，Excel中隐含有一个函数NumberString，它的功能是将数字转为中文大写数字，当然，也没有直接的逆函数，将中文大写数字转换为小写数字。
很想，将Word中的域代码:`{ SET myBK "一九一一年一月{ STYLEREF 1 \s }日" }{ myBK \@ "D" }`命名为`StringNumber`函数。

三、题注

Word默认插入带标题样式的题注，其域代码为图 `{ STYLEREF 1 \s }—{ SEQ 图 \* ARABIC \s 1 }`，很显然，我们只要将其中的`{ STYLEREF 1 \s }`替换为`{ SET myBK "一九一一年一月{ STYLEREF 1 \s }日" }{ myBK \@ "D" }`即可在中文大写数字的章标题中实现题注引用的常规数字（阿拉伯数字）。

替换方法：

1. 将域代码`{ SET myBK "一九一一年一月{ STYLEREF 1 \s }日" }{ myBK \@ "D" }`（即我命名的`StringNumber`函数）复制到剪贴板中；
2. 在正文中，按下ALT+F9，切换到域代码视图下（域的查找与替换，必须在显示域代码的情况下进行）。
3. 按下CTRL+H组合键，打开查找和替换对话框，在替换选项卡中，设置查找内容为“`^d STYLEREF 1 \s`”，在替换为中输入“^c”，不区分大小写，注意，“^d”后有一个半角空格，全部替换即可。

以下是引用sylun在2008-4-8 23:57:40的发言：

初步还有一种感觉，直接用Quote域而不用变量书签好像也行。
谢谢sylun兄的提醒！
用QUOTE域也行，并且更简单一些。QUOTE域对于单元格地址的引用也行（我以前的贴子），我想对于日期的引用同样应该没有问题！
QUOTE域是个宝！

简化后的域代码为：

注意,"{}"是由Ctrl+F9组合键自动插入的域标志!

域代码:

* 图 `{ QUOTE "一九一一年一月{ STYLEREF 1 \s }日" \@"D" }—{ SEQ 图 \* ARABIC \s 1 }`
* 表 `{ QUOTE "一九一一年一月{ STYLEREF 1 \s }日" \@"D" }—{ SEQ 表 \* ARABIC \s 1 }`
* 公式 `{ QUOTE "一九一一年一月{ STYLEREF 1 \s }日" \@"D" }—{ SEQ 公式 \* ARABIC \s 1 }`

借楼上几位高手的思路,做了一个宏,把下面的宏加到代码里,点功能区的insert caption按钮,就直接生成 `图1-1`这种形式的了.注意宏名字不能修改.

Word 2007英文版下测试通过,其他环境没测试.

```
Sub InsertCaption()
Dim ZH1 As String, ZH2 As String
Dim TH As String

Selection.TypeText "图"

ZH1 = "QUOTE ""一九一一年一月日"" \@ ""D"""
Selection.Fields.Add Range:=Selection.Range, PreserveFormatting:=False, Text:=ZH1
Selection.EndKey

ActiveWindow.View.ShowFieldCodes = True
Selection.MoveLeft , 11

ZH2 = "STYLEREF 1 \s"
Selection.Fields.Add Range:=Selection.Range, PreserveFormatting:=False, Text:=ZH2
Selection.EndKey
        
Selection.TypeText "-"

TH = "SEQ 图 \* ARABIC \s 1"
Selection.Fields.Add Range:=Selection.Range, PreserveFormatting:=False, Text:=TH

ActiveWindow.View.ShowFieldCodes = False
Selection.WholeStory
Selection.Fields.Update

Selection.EndKey

End Sub
```

## 页眉自动根据章节标题自动插入
```
{STYLEREF "标题 1"\n  \* MERGEFORMAT} {STYLEREF  "标题 1"  \* MERGEFORMAT}
```

## See also
* http://office-qa.com/Word/wd89.htm
* http://www.addbalance.com/usersguide/fields.htm
