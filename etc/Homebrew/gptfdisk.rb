class Gptfdisk < Formula
  desc "Text-mode GPT partitioning tools"
  homepage "https://sourceforge.net/projects/gptfdisk/"
  url "https://downloads.sourceforge.net/project/gptfdisk/gptfdisk/1.0.3/gptfdisk-1.0.3.tar.gz"
  sha256 "89fd5aec35c409d610a36cb49c65b442058565ed84042f767bba614b8fc91b5c"

  option "with-icu4c", "Use icu4c instead of internal functions for UTF-16 support. Use this if you are having problems with the new UTF-16 support."
  option "with-sgdisk", "Compile sgdisk."
  option "without-cgdisk", "Do not compile cgdisk."
  option "without-fixparts", "Do not compile fixparts."

  depends_on "icu4c" => :optional
  depends_on "popt" if build.with?("sgdisk")

  def install
    # Patch, upstream looks for wrong ncurses library
    inreplace "Makefile.mac", "/opt/local/lib/libncurses.a", "/usr/lib/libncurses.dylib"

    # Optional UTF-16 support from icu4c
    if build.with? "icu4c"
      inreplace "Makefile.mac", "-Wall", "-Wall -D USE_UTF16"
    end

    opts = ["gdisk"]
    opts << "sgdisk" if build.with? "sgdisk"
    opts << "cgdisk" if build.with? "cgdisk"
    opts << "fixparts" if build.with? "fixparts"

    system "make", "-f", "Makefile.mac", *opts
    sbin.install "gdisk"
    sbin.install "cgdisk" if build.with? "cgdisk"
    sbin.install "sgdisk" if build.with? "sgdisk"
    sbin.install "fixparts" if build.with? "fixparts"

    man8.install Dir["*.8"]
    doc.install Dir["*.html"]
  end

  test do
    assert_match /GPT fdisk \(gdisk\) version #{Regexp.escape(version)}/,
                 pipe_output("#{sbin}/gdisk", "\n")
  end
end
