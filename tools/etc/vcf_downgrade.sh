# If you are trying to view VCF 4.2 files in IGV - you may run into issues. This function might help you.
# This script will:
# 1. Rename the file as version 4.1
# 2. Replace parentheses in the INFO lines (IGV doesn't like these!)

function vcf_downgrade() {
  outfile=${1/.bcf/}
  outfile=${outfile/.gz/}
  outfile=${outfile/.vcf/}
  bcftools view --max-alleles 2 -O v $1 | \
  sed "s/##fileformat=VCFv4.2/##fileformat=VCFv4.1/" | \
  sed "s/(//" | \
  sed "s/)//" | \
  sed "s/,Version=\"3\">/>/" | \
  bcftools view -O z > ${outfile}.dg.vcf.gz
  tabix ${outfile}.dg.vcf.gz
}

