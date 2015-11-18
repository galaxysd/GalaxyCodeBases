# https://gist.github.com/goodbest/c39b0deef267552e2069
# https://twitter.com/lovegoodbest/status/666966149653995521
# 我再也受不了Office 2016 Mac 这种把相同的字体重复放在不同的程序包里的做法了...
# 于是把Excel、PowerPoint、OneNote、Outlook的字体文件夹直接软链到Word里的了。

# 补充：前两者的font和word的完全相同；后两者的font是word的子集。
# 运行之后没问题的话就可以把~/font_backup_xxxx这几个目录删了。
# 每次升级补丁之后应该需要重新运行一遍。

#Use word.app font path as real path
#"excel, powerpoint, onenote, outlook" font path are soft-linked to word.app's font path
basePATH="/Applications/"
WordPATH="Microsoft Word.app"
ExcelPATH="Microsoft Excel.app"
OutlookPATH="Microsoft Outlook.app"
PowerPointPATH="Microsoft PowerPoint.app"
OneNotePATH="Microsoft OneNote.app"
fontPATH="/Contents/Resources/Fonts"

#do Excel
sudo mv "$basePATH$ExcelPATH$fontPATH" ~/font_backup_Excel
sudo ln -s "$basePATH$WordPATH$fontPATH" "$basePATH$ExcelPATH$fontPATH"

#do Outlook
sudo mv "$basePATH$OutlookPATH$fontPATH" ~/font_backup_Outlook
sudo ln -s "$basePATH$WordPATH$fontPATH" "$basePATH$OutlookPATH$fontPATH"

#do PowerPoint
sudo mv "$basePATH$PowerPointPATH$fontPATH" ~/font_backup_PowerPoint
sudo ln -s "$basePATH$WordPATH$fontPATH" "$basePATH$PowerPointPATH$fontPATH"

#do OneNote
sudo mv "$basePATH$OneNotePATH$fontPATH" ~/font_backup_OneNote
sudo ln -s "$basePATH$WordPATH$fontPATH" "$basePATH$OneNotePATH$fontPATH"


echo "If everything is OK, you can delete folder: ~/font_bakcup_xxxx"