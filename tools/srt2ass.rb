#!/usr/bin/ruby
# == Synopsis
#
# str_ass: converts srt files in directory into ass files
#
# == Usage
#
# srt_ass.rb [OPTIONS]
#
# --help, -h
#    show help.
#
# --directory .(-d .)
#    directory to parse. current is default
#
# --mask anime_#(-m anime_#)
#    a mask which will be used in files replacement
require "getoptlong"
require 'rdoc/usage'
require "pp"
$KCODE = "u"

opts = GetoptLong.new(
  [ '--help', '-h', GetoptLong::NO_ARGUMENT ],
  [ '--directory', '-d', GetoptLong::OPTIONAL_ARGUMENT ],
  [ '--mask', '-m', GetoptLong::OPTIONAL_ARGUMENT ]
)

dir = "."
mask = nil

opts.each do |opt, arg|
  case opt
    when '--help'
      RDoc::usage

    when '--directory'
      dir = arg.to_s

    when '--mask'
      mask = arg.to_s.sub("#", '\\\1')
  end
end

# ass data
ass_header = "[Script Info]
Title: <untitled>
Original Script: <unknown>
Script Type: v4.00+
PlayResX: 0
PlayResY: 0
PlayDepth: 0

[V4+ Styles]
Format: Name, Fontname, Fontsize, PrimaryColour, SecondaryColour, OutlineColour, BackColour, Bold, Italic, Underline, StrikeOut, ScaleX, ScaleY, Spacing, Angle, BorderStyle, Outline, Shadow, Alignment, MarginL, MarginR, MarginV, Encoding
Style: Default,Flibustier,22,&H00FFFFFF,&H000000FF,&H0016360E,&H0017460B,0,0,0,0,100,100,0,0,1,1.4,0.6,2,10,10,10,1

[Events]
Format: Layer, Start, End, Style, Actor, MarginL, MarginR, MarginV, Effect, Text
"

#Style: Default,Comic Sans MS,23,&H00C3FFD2,&H0000FFFF,&H00000000,&H00000000,-1,0,0,0,100,100,0,0,1,3,1,2,10,10,12,204


ass_line = 'Dialogue: 0,\1.\2,\3.\4,Default,,0000,0000,0000,,\5'+"{ENDOFTEXT}\n"+'\6'
srt_line = /^\d+\n(\d\d:\d\d:\d\d),(\d\d)\d[\ ]?-->[\ ]?(\d\d:\d\d:\d\d),(\d\d)\d(?:  SSA.*)?\n([\s\S]+?)\n(^\d+\n\d\d:\d\d:\d\d)/


# get files
Dir.foreach(dir) do |f|
  next if [".", ".."].include?(f) || !f.match(/\.srt$/)

  data = nil
  File.open(dir+'/'+f, 'r') {|file| data = file.read }
  data = data.gsub("\r\n", "\n").gsub(/[\s\n]+$/m, '')+"\n0\n00:00:00"

  while data.match(srt_line) do
    data.gsub(srt_line) do |entry|
      tmp = entry.sub(srt_line, ass_line).split("{ENDOFTEXT}")
      tmp[0] = tmp[0].gsub("\n", "\\N").gsub(/\\N$/, '')
      new_entry = tmp.join("")

      data = data.sub(entry, new_entry)
      break
    end
  end
  if mask
    new_filename = (dir == "." ? "" : dir)+f.sub(/\.srt$/ ,'').sub(/.*\b(\d{4}|\d{3}|\d{2})\b.*/, mask)
    new_filename = new_filename+".ass" unless new_filename.match(/\.ass$/)
  else
    new_filename = dir+'/'+f.sub(/\.srt$/, '.ass')
  end
  p new_filename
  File.open(new_filename, 'w') {|file| file << ass_header << data.gsub(/\n0\n00:00:00$/m, '') }
end
