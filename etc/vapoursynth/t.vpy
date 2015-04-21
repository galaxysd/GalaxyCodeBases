import vapoursynth as vs
core = vs.get_core()

# only needed when libf3kdb.dylib wasn't installed
#core.std.LoadPlugin("/usr/local/Cellar/f3kdb/HEAD/lib/libf3kdb.dylib")
#core.std.LoadPlugin("/Users/Galaxy/git/ffms2/libffms2.4.dylib")

# http://forum.doom9.org/archive/index.php/t-165771-p-4.html
def fadeout(theclip):
    fade_input = theclip
    fadetime = 0.5
    fadeframes = fadetime * theclip.fps_num / theclip.fps_den
    beginframes = theclip.num_frames-1 - fadeframes
    blank = core.std.BlankClip(clip=fade_input, length=1, color=[0, 128, 128]) # assumes you're working with 8 bit YUV
    fade_frames = []
    for i in range(fade_input.num_frames):
        if (i <= beginframes) :
            fade_frames.append(fade_input[i])
        else :
            fade_frames.append(core.std.Merge(clipa=fade_input[i], clipb=blank, weight=[(i-beginframes)/fadeframes]))
    fade_output = core.std.Splice(clips=fade_frames)
    return fade_output

def fadein(theclip):
    fade_input = theclip
    fadetime = 1
    fadeframes = fadetime * theclip.fps_num / theclip.fps_den
    blank = core.std.BlankClip(clip=fade_input, length=1, color=[0, 128, 128]) # assumes you're working with 8 bit YUV
    fade_frames = []
    for i in range(fade_input.num_frames):
        if (i > fadeframes) :
            fade_frames.append(fade_input[i])
        else :
            fade_frames.append(core.std.Merge(clipa=fade_input[i], clipb=blank, weight=[(fadeframes-i)/fadeframes]))
    fade_output = core.std.Splice(clips=fade_frames)
    return fade_output

p1 = core.ffms2.Source(source='2.mp4')
p1 = core.std.Turn180(p1)
#p1 = core.text.ClipInfo(p1)
p1 = core.std.AddBorders(p1,top=124,bottom=124) # 568x320
#p1 = fadein(p1)
p1 = fadeout(p1)
#p1 = core.text.ClipInfo(p1)

p2 = core.ffms2.Source(source='3.mp4')
p2 = core.std.Turn180(p2)
#p2 = core.text.ClipInfo(p2)
p2 = core.std.AddBorders(p2,top=124,bottom=124)
p2 = fadein(p2)
p2 = fadeout(p2)

p3 = core.ffms2.Source(source='9.mp4')
p3 = core.std.Transpose(p3)
p3 = core.std.FlipHorizontal(p3)
#p3 = core.text.ClipInfo(p3)
p3 = core.std.AddBorders(p3,left=124,right=124)
p3 = fadein(p3)
#p3 = fadeout(p3)


#ret = core.std.Transpose(ret)
ret = core.std.Splice([p1,p2,p3])
ret.set_output()

# find . -name '*.mp4'|while read a;do echo "file $a"; done > mylist.txt
# ffmpeg -f concat -i mylist.txt -c copy output.mp4
# vspipe --y4m t.vpy - | mpv -
# vspipe --y4m t.vpy - | x264 --demuxer y4m - --output encoded.mkv
# ffmpeg -i output.mp4 -i encoded.mkv # see the list of stream numbers
# ffmpeg -i output.mp4 -i encoded.mkv -map 1:0 -map 0:1 -acodec copy -vcodec copy -shortest mg.mp4
# The -shortest option will cause the output duration to match the duration of the shortest input stream. http://ffmpeg.org/ffmpeg.html#Advanced-options

# proxychains4 -q brew tap mpv-player/mpv
# brew install --with-vapoursynth mpv
# proxychains4 git clone https://github.com/FFMS/ffms2.git
# ./configure --enable-shared && make
# cp ./src/core/.libs/libffms2.4.dylib /usr/local/Cellar/vapoursynth/26/lib/vapoursynth/
