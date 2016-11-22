class Rar < Formula
  desc "archive manager for RAR/ZIP file formats"
  homepage "http://www.rarlab.com"
  url "http://www.rarlab.com/rar/rarosx-5.4.0.tar.gz"
  version "5.40"
  sha256 "09a14f40718c68fc1c24b30acb55d0f2f90f3e13b372c48b6ef1e789d748b754"
  version_scheme 1

  bottle :unneeded

  resource "man" do
    url "http://manpages.ubuntu.com/manpages.gz/trusty/man1/rar.1.gz"
    sha256 "9a2dd38c3ec1e098f29b720116ad77a71e700dccebaa936cb7ace397702b287a"
  end

  def install
    bin.install "rar"
    lib.install "default.sfx"
    etc.install "rarfiles.lst"
    doc.install "acknow.txt", "order.htm", "rar.txt", "whatsnew.txt"
    man1.install resource("man")
  end

  test do
    cp test_fixtures("test.wav"), "test_orig.wav"
    system bin/"rar", "a", "test.rar", "test_orig.wav"
    system bin/"rar", "rn", "test.rar", "test_orig.wav", "test.wav"
    assert_match "test.wav", shell_output("#{bin}/rar l test.rar")
    system bin/"rar", "x", "test.rar"
    cmp "test.wav", "test_orig.wav"
  end
end
