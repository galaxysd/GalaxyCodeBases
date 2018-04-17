class Samtools01 < Formula
#class SamtoolsAT01 < Formula
  desc "Tools for manipulating next-generation sequencing data"
  homepage "https://samtools.sourceforge.io/"
  # doi "10.1093/bioinformatics/btp352"
  # tag "bioinformatics"
  # <http://samtools.sourceforge.net/pileup.shtml>
  # The `pileup' command has been removed in [0.1.17](https://github.com/samtools/samtools/commits/0.1.17) in commit c350a570e955d5b4c2c13f607b7442df6d332c67.
  # See <https://sourceforge.net/projects/samtools/files/samtools/0.1.16/> for README.

  url "https://github.com/samtools/samtools/archive/0.1.16.tar.gz"
  sha256 "7657e5dc66fbd1f02133349de34955a50635a9980c21ebff06c116e4e1e65986"

  keg_only :versioned_formula

  option "with-dwgsim", "Build with 'Whole Genome Simulation'"
  option "without-bcftools", "Do not install BCFtools"

  unless OS.mac?
    depends_on "ncurses"
    depends_on "zlib"
  end

  resource "dwgsim" do
    # http://sourceforge.net/apps/mediawiki/dnaa/index.php?title=Whole_Genome_Simulation
    url "https://downloads.sourceforge.net/project/dnaa/dwgsim/dwgsim-0.1.11.tar.gz"
    sha256 "6ffc8a4f7d20bc7c8b3efa1d2b3ae6cbf9609a93db976d4e7ccd2a209a2305b5"
  end

  stable do
    # https://github.com/lh3/samtools/issues/15
    # diff -uN ksort.h.0 ksort.h >> ../samtools-0.1.rb
    patch :DATA
  end

  def install
    system "make"
    system "make", "-C", "bcftools" if build.with? "bcftools"

    if build.with? "dwgsim"
      ohai "Building dwgsim"
      resource("dwgsim").stage do
        ln_s buildpath, "samtools"
        system "make", "CC=#{ENV.cc}"
        bin.install %w[dwgsim dwgsim_eval]
      end
    end

    bin.install %w[
      samtools bcftools/bcftools bcftools/vcfutils.pl
      misc/maq2sam-long misc/maq2sam-short misc/md5fa misc/md5sum-lite misc/wgsim
    ]

    bin.install Dir["misc/*.pl"]
    lib.install "libbam.a"
    man1.install "samtools.1"
    (share+"samtools").install "examples"
    (include+"bam").install Dir["*.h"]
  end

  test do
    assert_match "samtools", shell_output("#{bin}/samtools 2>&1", 1)
  end
end

__END__
diff --git a/ksort.h b/ksort.h
index fa850ab..f8d8c4c 100644
--- a/ksort.h
+++ b/ksort.h
@@ -141,7 +141,7 @@
 			tmp = *l; *l = l[i]; l[i] = tmp; ks_heapadjust_##name(0, i, l); \
 		}																\
 	}																	\
-	inline void __ks_insertsort_##name(type_t *s, type_t *t)			\
+	static inline void __ks_insertsort_##name(type_t *s, type_t *t)			\
 	{																	\
 		type_t *i, *j, swap_tmp;										\
 		for (i = s + 1; i < t; ++i)										\
