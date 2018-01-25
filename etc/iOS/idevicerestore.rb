class Idevicerestore < Formula
  desc "Tool that allows to restore firmware files to iOS devices."
  homepage "https://github.com/libimobiledevice/idevicerestore"
  url "https://github.com/libimobiledevice/idevicerestore/archive/master.tar.gz"
  version "0.0.1"
  #sha256 "dc174ad0bd1ddc498d697fa0298d1a97fecf50c93738af3dadcb31e74c6f5748"

  depends_on "automake" => :build
  depends_on "libtool" => :build
  depends_on "autoconf" => :build
  depends_on "pkg-config" => :build
  depends_on "libirecovery"
  depends_on "libimobiledevice"
  depends_on "libzip"

  def install
    system "./autogen.sh", "--disable-debug",
                          "--disable-dependency-tracking",
                          "--disable-silent-rules",
                          "--prefix=#{prefix}"
    system "make"
    system "make", "install"
  end

  test do
    system "#{bin}/idevicerestore", "--help"
  end
end
