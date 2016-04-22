IndivID: An Application for Identification of Individuals with Genotypes from Multiple Microsatellite Loci

As microsatellites often present high levels of polymorphisms within species, they are often used as molecular markers for intra-specific study. In most genetic-based studies of wildlife, non-invasive sampling of genetic materials, such as hair, feces, feathers, and so on, is extensively used to get animal DNA in the research field of population genetics, molecular ecology, conservation biology, and etc. Because researchers don't handle animals directly, one individual is always sampled repeatedly. So identification of individuals is necessary for performing subsequent analysis. Several loci of microsatellite are often informative enough to identify individuals and get population genetic profile. IndivID is a program designed for identification of individuals using samples with multiple microsatellite loci. The program is written in Perl.

1. Criterion
If genotyping of a locus is repeated three times, false allele hardly occurs, but allele dropout can’t be avoided. So we define one locus of two samples as “Compatible Locus (CL)” if there is no false allele between them. For example, if the diploid genotypes of two samples of a locus are “a/a, a/b” or “a/a, b/b”, the locus of the two samples are CL; if one or two of the two samples are missing data, the locus are also CL; but if more than two alleles appear in one locus of two samples, like “a/a, b/c” or “a/b, b/c” or “a/b, c/d”, the locus are not CL. Then we define two samples as “Compatible Samples (CS)” if all loci of the two samples are CL. After defining every two samples as CS or not CS, we infer that a group of samples belong to one individual, if every two samples in the group are CS and every sample in the group with any sample out of the group are not CS; we define a group of samples as uncertain groups, if not every two samples in the group are CS and every sample in the group with any sample out of the group are not CS. Besides the inferred individuals and the defined uncertain groups, we infer that each of the remained samples is one individual.

2. Input
The program can only deal with diploid genotype now. The input data for each sample must be stored in one row, where each locus is in two consecutive columns. And the columns must be tab-delimited. There must be one pre-genotype column which records sample ID and one pre-genotype row which records locus ID. The loci IDs are delimited by double tabs.

3. Output
The program outputs one file which contains three part. The first part is inferred individuals. For each individual, the first row is inferred genotype, and the following rows are genotypes of samples belong to the individual. The second part is uncertain groups. Samples belong to each group are listed below the group ID. The third part contains individuals with only one sample.

4. Run the program
The program is written in Perl, so the operating system must have Perl installed. Copy the program code file and the input data file to an appropriate folder. Then open a command-line interface. Run the program by typing “perl <code file name> <input file name> <output file name>” and press enter.

5. Availability
IndivID is written in Perl 5. The source code and documentation can be downloaded freely at: https://github.com/WoodyMiao/IndivID
