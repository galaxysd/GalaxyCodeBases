class Idevicelocation < Formula
  desc "Library to communicate with iOS devices natively"
  homepage "https://github.com/JonGabilondoAngulo/idevicelocation"
  url "https://github.com/galaxy001/idevicelocation/archive/v101.tar.gz"
  sha256 "cfaf45a431d872676c688389bd25a6399da210c243971203b32e5bcb7de6b641"
  license "LGPL-2.1"

  depends_on "autoconf"
  depends_on "automake"
  depends_on "libtool"
  depends_on "pkg-config"
  depends_on "libplist"
  depends_on "libimobiledevice"
  depends_on "libzip"
  depends_on "openssl@1.1"

  def install
    system "./autogen.sh", "--disable-dependency-tracking",
                           "--disable-silent-rules",
                           "--prefix=#{prefix}"
    system "make", "install"
  end

  test do
    system "#{bin}/idevicelocation", "--help"
  end
end
