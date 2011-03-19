package Galaxy::Data;
#package main;
use strict;
require Exporter;
our @ISA   =qw(Exporter);
our @EXPORT    =qw(%IUB %REV_IUB %IUP %gen_code);
our @EXPORT_OK   =qw(%Alphabets %Alphabets_strict);
our $VERSION   = v1.0.0;

=head1 Purpose
Provide some common data.
=cut

#http://doc.bioperl.org/bioperl-live/Bio/Tools/IUPAC.html#BEGIN1
our %IUB = ( A => [qw(A)],
	     C => [qw(C)],
	     G => [qw(G)],
	     T => [qw(T)],
	     U => [qw(U)],
	     M => [qw(A C)],
	     R => [qw(A G)],
	     W => [qw(A T)],
	     S => [qw(C G)],
	     Y => [qw(C T)],
	     K => [qw(G T)],
	     V => [qw(A C G)],
	     H => [qw(A C T)],
	     D => [qw(A G T)],
	     B => [qw(C G T)],
	     X => [qw(G A T C)],
	     N => [qw(G A T C)]
	     );
our %REV_IUB = (A	=> 'A',
		T	=> 'T',
		C	=> 'C',
		G 	=> 'G',
		AC	=> 'M',
		AG	=> 'R',
		AT	=> 'W',
		CG	=> 'S',
		CT	=> 'Y',
		'GT'	=> 'K',
		ACG	=> 'V',
		ACT	=> 'H',
		AGT	=> 'D',
		CGT	=> 'B',
		ACGT=> 'N',
		N	=> 'N'
		);

our %IUP = (A => [qw(A)],
	    B => [qw(D N)],
	    C => [qw(C)],
	    D => [qw(D)],
	    E => [qw(E)],
	    F => [qw(F)],
	    G => [qw(G)],
	    H => [qw(H)],
	    I => [qw(I)],
        J => [qw(I L)],
	    K => [qw(K)],
	    L => [qw(L)],
	    M => [qw(M)],
	    N => [qw(N)],
        O => [qw(O)],
	    P => [qw(P)],
	    Q => [qw(Q)],
	    R => [qw(R)],
	    S => [qw(S)],
	    T => [qw(T)],
	    U => [qw(U)],
	    V => [qw(V)],
	    W => [qw(W)],
	    X => [qw(X)],
	    Y => [qw(Y)],
	    Z => [qw(E Q)],
	    '*' => ['*']
		);

#http://doc.bioperl.org/releases/bioperl-1.4/Bio/Tools/SeqStats.html#BEGIN1
our %Alphabets =   (
		    'dna'     => [ qw(A C G T R Y M K S W H B V D X N) ],
		    'rna'     => [ qw(A C G U R Y M K S W H B V D X N) ],
		    'protein' => [ qw(A R N D C Q E G H I L K M F U
				      P S T W X Y V B Z *) ], # sac: added B, Z
		    );

# SAC: new strict alphabet: doesn't allow any ambiguity characters.
our %Alphabets_strict = (
			 'dna'     => [ qw( A C G T ) ],
			 'rna'     => [ qw( A C G U ) ],
			 'protein'    => [ qw(A R N D C Q E G H I L K M F U
					      P S T W Y V) ],
			 );
#  IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE:
#   Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030.

# From GuoXS
our %gen_code = ( "TTT" => "Phe", "TTC" => "Phe", "TTA" => "Leu", "TTG" => "Leu", "CTT" => "Leu", "CTC" => "Leu", "CTA" => "Leu", "CTG" => "Leu", "ATT" => "Ile", "ATC" => "Ile", "ATA" => "Ile", "ATG" => "Met", "GTT" => "Val", "GTC" => "Val", "GTA" => "Val", "GTG" => "Val", "TCT" => "Ser", "TCC" => "Ser", "TCA" => "Ser", "TCG" => "Ser", "CCT" => "Pro", "CCC" => "Pro", "CCA" => "Pro", "CCG" => "Pro", "ACT" => "Thr", "ACC" => "Thr", "ACA" => "Thr", "ACG" => "Thr", "GCT" => "Ala", "GCC" => "Ala", "GCA" => "Ala", "GCG" => "Ala",  "TAT" => "Tyr", "TAC" => "Tyr", "TAA" => "STOP", "TAG" => "STOP", "CAT" => "His", "CAC" => "His", "CAA" => "Gln", "CAG" => "Gln", "AAT" => "Asn", "AAC" => "Asn", "AAA" => "Lys", "AAG" => "Lys", "GAT" => "Asp", "GAC" => "Asp", "GAA" => "Glu", "GAG" => "Glu",  "TGT" => "Cys", "TGC" => "Cys", "TGA" => "STOP", "TGG" => "Trp", "CGT" => "Arg", "CGC" => "Arg", "CGA" => "Arg", "CGG" => "Arg", "AGT" => "Ser", "AGC" => "Ser", "AGA" => "Arg", "AGG" => "Arg", "GGT" => "Gly", "GGC" => "Gly", "GGA" => "Gly", "GGG" => "Gly" );

1;
