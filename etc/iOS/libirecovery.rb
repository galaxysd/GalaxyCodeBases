class Libirecovery < Formula
  desc "library for communication to iBoot/iBSS on iOS devices via USB."
  homepage "https://github.com/libimobiledevice/libirecovery"
  url "https://github.com/libimobiledevice/libirecovery/archive/master.tar.gz"
  version "0.2.0"
  #sha256 "f5b685aced24ad3ea4f06a378e0f4582a69a4d39750f6e529ba2fff212785ccc"

  depends_on "automake" => :build
  depends_on "libtool" => :build
  depends_on "autoconf" => :build

  def install
    system "./autogen.sh", "--disable-debug",
                          "--disable-dependency-tracking",
                          "--disable-silent-rules",
                          "--prefix=#{prefix}"
    system "make"
    system "make", "install"
  end

  test do
    system "ls", "#{lib}/libirecovery.dylib"
  end
end
