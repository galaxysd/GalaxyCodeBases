http://sourceforge.net/projects/anytag/files/anytag2.0/

anytag-2.0

[2.0 Main Changes]

1, Optimize the short reads alignment.
1.1, Remove huge temporary files (which store the intermedia hits)
1.2, Load all index in memory, so that one read can query all indexs at once
1.3, Number of seeds can be customized by user

2, Optimize local assembly
2.1, Add re-alignment, retrieve more reads in local group. Thus, increase the accuracy of consensus sequence
2.2, Remove the pair_score module, I plan to estimate pair_score in future
2.3, Add MSA and SAM file for all FIS. MSA file is multiple alignments.
2.4, Filter read pairs which only have one end inside FIS
2.5, Remove PCR duplication

3, Add new command `all`, to more efficently invoker aln and asm

4, Add new function `lnk`. I try to use mate-paired reads to build mate information for FIS. It is in development.

[The purpose of anytag]

NGS can product a read pairs which come from one original DNA fragment. If the internal gap can be correctly filled,
we will get the full length sequence of original DNA fragment.

[FIS]

FIS is short for Filled-In Sequence. The sequence of a FIS repreasents the original DNA.

[AR]

AR is short for Anchoring Reads.  If the gap between two ends of one AR is well filled, it become FIS.

[SR]

SR is short for Supporting Reads. SRs are used to fill in the gap of AR. The insert-size of SRs should be less than AR.

For example, all of our reads are 2*100bp,
1, inserts(bp): [250(ARs)], there is only ARs, the length of FIS will be 250bp
2, inserts(bp): [200(SRs), 400(ARs)], the length of FIS will be 400bp
3, inserts(bp): [200(SRs), 300(SRs), 600(ARs)], the length of FIS will be 600bp
4, inserts(bp): [200, 250, 300, 350, 400, 450, 500, 550(|<-SRs),600(ARs)], the length of FIS will be 600bp
5, inserts(bp): [200, 250, 300, 350, 400, 450, 500, 550, 600 (|<-SRs), 1000(ARs)], requiring much more sequence depth
                  the length of FIS will be 1000bp

[How does anytag work]

1, Clustering. Align SRs to AR. The minimum overlap between AR and SRs should be as large as to be nearly unique, such as 30bp.
   All alignments are forword. That is, AR and SRs are the same strand.
2, Local assembly. Build an overlap graph for AR and SRs. Pairwise alignment for all AR and SRs is done to find overlaps.
   Then, find a path from one end of AR to the other end. To make sure the path is correct, it prior to traverse max overlaps.
   Final, call consensus sequence (FIS) for AR.

[Can I use anytag in a huge genome]

Yes, anytag-2.0 can handle huge genome. It is fast enough. Roughly, constructing 10X FIS for human require 32 cpu * 3 days.

[Repeats]

If a repeat is less than the insert size of AR, anytag can build correct FISs which crossing it.

[Why FIS is better for whole genome assembly than direct short reads assembly]

1, for repeats. anytag can solve most of repeats of size less than the insert size of AR.
2, for heterozygosity. anytag use smith-waterman algorithm to align short reads in local assembly. Most of short reads assemblers cannot.

In our simulations of dm3 genome (50X, 2 * 80bp, max insert size 570bp):
------------------------------
|Heterozygosity|Method|N50   |
------------------------------
|0.001         |anytag|150.0k|
|0.001         |velvet|16.88k|
------------------------------
|0.01          |anytag|136.0k|
|0.01          |velvet|3.48k |
------------------------------
|0.02          |anytag|125.6k|
|0.02          |velvet|1.81k |
------------------------------

[The accuracy of FIS]

The accuracy of FIS is higher than first several bases of solexa reads, also higher than sanger reads.

[FIS in resequencing]

Although I haven't tested it in resequencing project, it will be useful in detecting small INDEL and SV, especially in detceting the breakpoints.

[Feedback]

ruanjue@gmail.com ruanj@big.ac.cn Beijing Institute of Genomics, Chinese Academy of Sciences

