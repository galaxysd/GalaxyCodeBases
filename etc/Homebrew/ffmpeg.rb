class Ffmpeg < Formula
  desc "Play, record, convert, and stream audio and video"
  homepage "https://ffmpeg.org/"
  url "https://ffmpeg.org/releases/ffmpeg-3.0.1.tar.bz2"
  sha256 "f7f7052c120f494dd501f96becff9b5a4ae10cfbde97bc2f1e9f0fd6613a4984"
  head "https://github.com/FFmpeg/FFmpeg.git"

  bottle do
    sha256 "754cd2da64d2dda3aab3bca329757e102c781678ee9289e3ed238f5a639229cf" => :el_capitan
    sha256 "b9a08881ea2053b739b92cb7dfdef4ce704ffef1776c68173e5e08c996c94004" => :yosemite
    sha256 "1ce859c211f6759e1eac392ca85f91f99d89e464520312edd9f96618e97626b3" => :mavericks
  end

  option "without-x264", "Disable H.264 encoder"
  option "without-lame", "Disable MP3 encoder"
  option "without-xvid", "Disable Xvid MPEG-4 video encoder"
  option "without-qtkit", "Disable deprecated QuickTime framework"

  option "with-rtmpdump", "Enable RTMP protocol"
  option "with-libass", "Enable ASS/SSA subtitle format"
  option "with-opencore-amr", "Enable Opencore AMR NR/WB audio format"
  option "with-openjpeg", "Enable JPEG 2000 image format"
  option "with-openssl", "Enable SSL support"
  option "with-libssh", "Enable SFTP protocol via libssh"
  option "with-schroedinger", "Enable Dirac video format"
  option "with-ffplay", "Enable FFplay media player"
  option "with-tools", "Enable additional FFmpeg tools"
  option "with-fdk-aac", "Enable the Fraunhofer FDK AAC library"
  option "with-libvidstab", "Enable vid.stab support for video stabilization"
  option "with-x265", "Enable x265 encoder"
  option "with-libsoxr", "Enable the soxr resample library"
  option "with-webp", "Enable using libwebp to encode WEBP images"
  option "with-zeromq", "Enable using libzeromq to receive commands sent through a libzeromq client"
  option "with-snappy", "Enable Snappy library"
  option "with-dcadec", "Enable dcadec library"
  option "with-rubberband", "Enable rubberband library"
  option "with-zimg", "Enable z.lib zimg library"
  option "with-openh264", "Enable OpenH264 library"

  depends_on "pkg-config" => :build

  # manpages won't be built without texi2html
  depends_on "texi2html" => :build
  depends_on "yasm" => :build

  depends_on "x264" => :recommended
  depends_on "lame" => :recommended
  depends_on "xvid" => :recommended

  depends_on "faac" => :optional
  depends_on "fontconfig" => :optional
  depends_on "freetype" => :optional
  depends_on "theora" => :optional
  depends_on "libvorbis" => :optional
  depends_on "libvpx" => :optional
  depends_on "rtmpdump" => :optional
  depends_on "opencore-amr" => :optional
  depends_on "libass" => :optional
  depends_on "openjpeg" => :optional
  depends_on "sdl" if build.with? "ffplay"
  depends_on "snappy" => :optional
  depends_on "speex" => :optional
  depends_on "schroedinger" => :optional
  depends_on "fdk-aac" => :optional
  depends_on "opus" => :optional
  depends_on "frei0r" => :optional
  depends_on "libcaca" => :optional
  depends_on "libbluray" => :optional
  depends_on "libsoxr" => :optional
  depends_on "libvidstab" => :optional
  depends_on "x265" => :optional
  depends_on "openssl" => :optional
  depends_on "libssh" => :optional
  depends_on "webp" => :optional
  depends_on "zeromq" => :optional
  depends_on "libbs2b" => :optional
  depends_on "dcadec" => :optional
  depends_on "rubberband" => :optional
  depends_on "zimg" => :optional
  depends_on "openh264" => :optional

  patch :DATA

  def install
    args = ["--prefix=#{prefix}",
            "--enable-shared",
            "--enable-pthreads",
            "--enable-gpl",
            "--enable-version3",
            "--enable-hardcoded-tables",
            "--enable-avresample",
            "--cc=#{ENV.cc}",
            "--host-cflags=#{ENV.cflags}",
            "--host-ldflags=#{ENV.ldflags}",
           ]

    args << "--enable-opencl" if MacOS.version > :lion

    args << "--enable-libx264" if build.with? "x264"
    args << "--enable-libmp3lame" if build.with? "lame"
    args << "--enable-libxvid" if build.with? "xvid"
    args << "--enable-libsnappy" if build.with? "snappy"

    args << "--enable-libfontconfig" if build.with? "fontconfig"
    args << "--enable-libfreetype" if build.with? "freetype"
    args << "--enable-libtheora" if build.with? "theora"
    args << "--enable-libvorbis" if build.with? "libvorbis"
    args << "--enable-libvpx" if build.with? "libvpx"
    args << "--enable-librtmp" if build.with? "rtmpdump"
    args << "--enable-libopencore-amrnb" << "--enable-libopencore-amrwb" if build.with? "opencore-amr"
    args << "--enable-libfaac" if build.with? "faac"
    args << "--enable-libass" if build.with? "libass"
    args << "--enable-ffplay" if build.with? "ffplay"
    args << "--enable-libssh" if build.with? "libssh"
    args << "--enable-libspeex" if build.with? "speex"
    args << "--enable-libschroedinger" if build.with? "schroedinger"
    args << "--enable-libfdk-aac" if build.with? "fdk-aac"
    args << "--enable-openssl" if build.with? "openssl"
    args << "--enable-libopus" if build.with? "opus"
    args << "--enable-frei0r" if build.with? "frei0r"
    args << "--enable-libcaca" if build.with? "libcaca"
    args << "--enable-libsoxr" if build.with? "libsoxr"
    args << "--enable-libvidstab" if build.with? "libvidstab"
    args << "--enable-libx265" if build.with? "x265"
    args << "--enable-libwebp" if build.with? "webp"
    args << "--enable-libzmq" if build.with? "zeromq"
    args << "--enable-libbs2b" if build.with? "libbs2b"
    args << "--enable-libdcadec" if build.with? "dcadec"
    args << "--enable-librubberband" if build.with? "rubberband"
    args << "--enable-libzimg" if build.with? "zimg"
    args << "--disable-indev=qtkit" if build.without? "qtkit"
    args << "--enable-libopenh264" if build.with? "openh264"

    if build.with? "openjpeg"
      args << "--enable-libopenjpeg"
      args << "--disable-decoder=jpeg2000"
      args << "--extra-cflags=" + `pkg-config --cflags libopenjpeg`.chomp
    end

    # These librares are GPL-incompatible, and require ffmpeg be built with
    # the "--enable-nonfree" flag, which produces unredistributable libraries
    if %w[faac fdk-aac openssl].any? { |f| build.with? f }
      args << "--enable-nonfree"
    end

    # A bug in a dispatch header on 10.10, included via CoreFoundation,
    # prevents GCC from building VDA support. GCC has no problems on
    # 10.9 and earlier.
    # See: https://github.com/Homebrew/homebrew/issues/33741
    if MacOS.version < :yosemite || ENV.compiler == :clang
      args << "--enable-vda"
    else
      args << "--disable-vda"
    end

    # For 32-bit compilation under gcc 4.2, see:
    # https://trac.macports.org/ticket/20938#comment:22
    ENV.append_to_cflags "-mdynamic-no-pic" if Hardware.is_32_bit? && Hardware::CPU.intel? && ENV.compiler == :clang

    system "./configure", *args

    if MacOS.prefer_64_bit?
      inreplace "config.mak" do |s|
        shflags = s.get_make_var "SHFLAGS"
        if shflags.gsub!(" -Wl,-read_only_relocs,suppress", "")
          s.change_make_var! "SHFLAGS", shflags
        end
      end
    end

    system "make", "install"

    if build.with? "tools"
      system "make", "alltools"
      bin.install Dir["tools/*"].select { |f| File.executable? f }
    end
  end

  def caveats
    if build.without? "faac" then <<-EOS.undent
      The native FFmpeg AAC encoder has been stable since FFmpeg 3.0. If you
      were using libvo-aacenc or libaacplus, both of which have been dropped in
      FFmpeg 3.0, please consider switching to the native encoder (-c:a aac),
      fdk-aac (-c:a libfdk_aac, ffmpeg needs to be installed with the
      --with-fdk-aac option), or faac (-c:a libfaac, ffmpeg needs to be
      installed with the --with-faac option).

      See the announcement
      https://ffmpeg.org/index.html#removing_external_aac_encoders for details,
      and https://trac.ffmpeg.org/wiki/Encode/AAC on best practices of encoding
      AAC with FFmpeg.
      EOS
    end
  end

  test do
    # Create an example mp4 file
    system "#{bin}/ffmpeg", "-y", "-filter_complex",
        "testsrc=rate=1:duration=1", "#{testpath}/video.mp4"
    assert (testpath/"video.mp4").exist?
  end
end

__END__
diff --git a/libavfilter/vf_owdenoise.c b/libavfilter/vf_owdenoise.c
index 3a77f89..6becccc 100644
--- a/libavfilter/vf_owdenoise.c
+++ b/libavfilter/vf_owdenoise.c
@@ -222,7 +222,6 @@ static void filter(OWDenoiseContext *s,
 
 static int filter_frame(AVFilterLink *inlink, AVFrame *in)
 {
-    int direct = 0;
     AVFilterContext *ctx = inlink->dst;
     OWDenoiseContext *s = ctx->priv;
     AVFilterLink *outlink = ctx->outputs[0];
@@ -231,8 +230,14 @@ static int filter_frame(AVFilterLink *inlink, AVFrame *in)
     const int ch = AV_CEIL_RSHIFT(inlink->h, s->vsub);
 
     if (av_frame_is_writable(in)) {
-        direct = 1;
         out = in;
+
+        if (s->luma_strength > 0)
+            filter(s, out->data[0], out->linesize[0], in->data[0], in->linesize[0], inlink->w, inlink->h, s->luma_strength);
+        if (s->chroma_strength > 0) {
+            filter(s, out->data[1], out->linesize[1], in->data[1], in->linesize[1], cw,        ch,        s->chroma_strength);
+            filter(s, out->data[2], out->linesize[2], in->data[2], in->linesize[2], cw,        ch,        s->chroma_strength);
+        }
     } else {
         out = ff_get_video_buffer(outlink, outlink->w, outlink->h);
         if (!out) {
@@ -240,13 +245,20 @@ static int filter_frame(AVFilterLink *inlink, AVFrame *in)
             return AVERROR(ENOMEM);
         }
         av_frame_copy_props(out, in);
-    }
 
-    filter(s, out->data[0], out->linesize[0], in->data[0], in->linesize[0], inlink->w, inlink->h, s->luma_strength);
-    filter(s, out->data[1], out->linesize[1], in->data[1], in->linesize[1], cw,        ch,        s->chroma_strength);
-    filter(s, out->data[2], out->linesize[2], in->data[2], in->linesize[2], cw,        ch,        s->chroma_strength);
+        if (s->luma_strength > 0) {
+            filter(s, out->data[0], out->linesize[0], in->data[0], in->linesize[0], inlink->w, inlink->h, s->luma_strength);
+        } else {
+            av_image_copy_plane(out->data[0], out->linesize[0], in ->data[0], in ->linesize[0], inlink->w, inlink->h);
+        }
+        if (s->chroma_strength > 0) {
+            filter(s, out->data[1], out->linesize[1], in->data[1], in->linesize[1], cw, ch, s->chroma_strength);
+            filter(s, out->data[2], out->linesize[2], in->data[2], in->linesize[2], cw, ch, s->chroma_strength);
+        } else {
+            av_image_copy_plane(out->data[1], out->linesize[1], in ->data[1], in ->linesize[1], inlink->w, inlink->h);
+            av_image_copy_plane(out->data[2], out->linesize[2], in ->data[2], in ->linesize[2], inlink->w, inlink->h);
+        }
 
-    if (!direct) {
         if (in->data[3])
             av_image_copy_plane(out->data[3], out->linesize[3],
                                 in ->data[3], in ->linesize[3],
