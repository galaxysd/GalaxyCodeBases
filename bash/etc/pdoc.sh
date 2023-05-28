#!/bin/sh

inPath="$1"
#inputname=${infile%%.*}
#allext="${infile#*.}"

mkdir -p ./work
#rm ./work/*.png
echo "# ${inPath}  " > "${inPath}.md"

find "./$inPath" -type f -exec bash -c '
for item do
  item=${item#foo}
  filename=$(basename -- "$item")
  newname="./work/'${inPath}'.$filename"
  if [[ ! -e "$newname" ]]; then
  	convert "$item" -colors 128 -type Palette -set units PixelsPerInch -density 315 "$newname"
  fi
done
' bash {} +

find "./work" -type f -name "${inPath}.*" -print0 | sort -z | xargs -r0 -I{} echo '![]('{}')  ' >> "${inPath}.md"
pandoc -f markdown -t docx "${inPath}.md" -o "${inPath}.docx"
