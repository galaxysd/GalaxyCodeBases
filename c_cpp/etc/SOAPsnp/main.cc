#include "soap_snp.h"
#include <getopt.h>

using namespace std;

int usage() {
	cerr<<"SoapSNP"<<endl;
	cerr<<"Compulsory Parameters:"<<endl;
	cerr<<"-i <FILE> Input SORTED Soap Result"<<endl;
	cerr<<"-d <FILE> Reference Sequence in fasta format"<<endl;
	cerr<<"-o <FILE> Output consensus file"<<endl;
	cerr<<"Optional Parameters:(Default in [])"<<endl;
	cerr<<"-z <Char> ASCII chracter standing for quality==0 [@]"<<endl;
	cerr<<"-g <Double> Global Error Dependency Coefficient, 0.0(complete dependent)~1.0(complete independent)[0.9]"<<endl;
	cerr<<"-p <Double> PCR Error Dependency Coefficient, 0.0(complete dependent)~1.0(complete independent)[0.5]"<<endl;
	cerr<<"-r <Double> novel altHOM prior probability [0.0005]"<<endl;
	cerr<<"-e <Double> novel HET prior probability [0.0010]"<<endl;
	cerr<<"-t set transition/transversion ratio to 2:1 in prior probability"<<endl;
	cerr<<"-s <FILE> Pre-formated dbSNP information"<<endl;
	cerr<<"-2 specify this option will REFINE SNPs using dbSNPs information [Off]"<<endl;
	cerr<<"-a <Double> Validated HET prior, if no allele frequency known [0.1]"<<endl;
	cerr<<"-b <Double> Validated altHOM prior, if no allele frequency known[0.05]"<<endl;
	cerr<<"-j <Double> Unvalidated HET prior, if no allele frequency known [0.02]"<<endl;
	cerr<<"-k <Double> Unvalidated altHOM rate, if no allele frequency known[0.01]"<<endl;
	cerr<<"-u Enable rank sum test to give HET further penalty for better accuracy. [Off]"<<endl;
	//cerr<<"-n Enable binomial probability calculation to give HET for better accuracy. [Off]"<<endl;
	cerr<<"-m Enable monoploid calling mode, this will ensure all consensus as HOM and you probably should SPECIFY higher altHOM rate. [Off]"<<endl;
	cerr<<"-q Only output potential SNPs. Useful in Text output mode. [Off]"<<endl;
	cerr<<"-M <FILE> Output the quality calibration matrix; the matrix can be reused with -I if you rerun the program"<<endl;
	cerr<<"-I <FILE> Input previous quality calibration matrix. It cannot be used simutaneously with -M"<<endl;
	cerr<<"-L <short> maximum length of read [45]"<<endl;
	cerr<<"-Q <short> maximum FASTQ quality score [40]"<<endl;
	cerr<<"-F <int> Output format. 0: Text; 1: GLFv2; 2: GPFv2.[0]"<<endl;
	cerr<<"-E <String> Extra headers EXCEPT CHROMOSOME FIELD specified in GLFv2 output. Format is \"TypeName1:DataName1:TypeName2:DataName2\"[""]"<<endl;
	cerr<<"-T <FILE> Only call consensus on regions specified in FILE. Format: ChrName\\tStart\\tEnd."<<endl;
	//cerr<<"-S <FILE> Output summary of consensus"<<endl;
	cerr<<"-h Display this help"<<endl;
	exit(1);
	return 0;
}

int readme() {
	return usage();
}

int main ( int argc, char * argv[]) {
	// This part is the default values of all parameters
	Parameter * para = new Parameter;
	std::string alignment_name, consensus_name("");
	bool is_matrix_in = false; // Generate the matrix or just read it?
	int c;
	Files files;
	while((c=getopt(argc,argv,"i:d:o:z:g:p:r:e:ts:2a:b:j:k:unmqM:I:L:Q:S:F:E:T:h")) != -1) {
		switch(c) {
			case 'i':
			{
				// Soap Alignment Result
				files.soap_result.clear();
				files.soap_result.open(optarg);
				if( ! files.soap_result) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					exit(1);
				}
				alignment_name = optarg;
				break;
			}
			case 'd':
			{
				// The reference genome in fasta format
				files.ref_seq.clear();
				files.ref_seq.open(optarg);
				if( ! files.ref_seq) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					exit(1);
				}
				files.ref_seq.clear();
				break;
			}
			case 'o':
			{
				consensus_name = optarg;
				break;
			}
			case 'z':
			{
				// The char stands for quality==0 in fastq format
				para->q_min = optarg[0];
				if(para->q_min == 33) {
					clog<<"Standard Fastq System Set"<<endl;
				}
				else if(para->q_min == 64) {
					clog<<"Illumina Fastq System Set"<<endl;
				}
				else {
					clog<<"Other types of Fastq files?? Are you sure?"<<endl;
				}
				para->q_max = para->q_min + 40;
				break;
			}
			case 'g':
			{
				para->global_dependency= log10(atof(optarg));
				break;
			}
			case 'p':
			{
				para->pcr_dependency= log10(atof(optarg));
				break;
			}
			case 'r':
			{
				para->althom_novel_r = atof(optarg);
				break;
			}
			case 'e':
			{
				para->het_novel_r=atof(optarg);
				break;
			}
			case 't':
			{
				para->transition_dominant=true;
				break;
			}
			case 's':
			{
				// Optional: A pre-formated dbSNP table
				files.dbsnp.clear();
				files.dbsnp.open(optarg);
				if( ! files.ref_seq) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					exit(1);
				}
				files.dbsnp.clear();
				break;
			}
			case '2':
			{
				// Refine prior probability based on dbSNP information
				para->refine_mode = true;
				break;
			}
			case 'a':
			{
				para->althom_val_r=atof(optarg);
				break;
			}
			case 'b':
			{
				para->het_val_r=atof(optarg);
				break;
			}
			case 'j':
			{
				para->althom_unval_r=atof(optarg);
				break;
			}
			case 'k':
			{
				para->het_unval_r=atof(optarg);
				break;
			}
			case 'u':
			{
				para->rank_sum_mode = true;
				break;
			}
			case 'n':
			{
				para->binom_mode = true;
				break;
			}
		 	case 'm':
			{
				para->is_monoploid=1;
				break;
			}
			case 'q':
			{
				para->is_snp_only=1;
				break;
			}
			case 'M':
			{
				files.matrix_file.close(); files.matrix_file.clear();
				// Output the calibration matrix
				files.matrix_file.open(optarg, fstream::out);
				if( ! files.matrix_file) {
					cerr<<"Cannot creat file :"<<optarg<<endl;
					exit(1);
				}
				files.matrix_file.clear();
				break;
			}
			case 'I':
			{
				files.matrix_file.close(); files.matrix_file.clear();
				// Input the calibration matrix
				files.matrix_file.open(optarg, fstream::in);
				if( ! files.matrix_file) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					exit(1);
				}
				files.matrix_file.clear();
				is_matrix_in = true;
				break;
			}
			case 'S':
			{
				//files.summary.open(optarg);
				//// Output the summary of consensus
				//if( ! files.summary ) {
				//	cerr<<"No such file or directory: "<<optarg<<endl;
				//	exit(1);
				//}
				break;
			}
			case 'L':
			{
				para->read_length = atoi(optarg);
				break;
			}
			case 'Q':
			{
				para->q_max = optarg[0];
				if(para->q_max < para->q_min) {
					cerr<< "FASTQ quality character error: Q_MAX > Q_MIN" <<endl;
				}
				break;
			}
			case 'F': {
				para->glf_format = atoi(optarg);
				break;
			}
			case 'E': {
				para->glf_header = optarg;
				break;
			}
			case 'T': {
				files.region.clear();
				files.region.open(optarg);
				files.region.clear();
				para->region_only = true;
				break;
			}
			case 'h':readme();break;
			case '?':usage();break;
			default: cerr<<"Unknown error in command line parameters"<<endl;
		}
	}
	if( consensus_name=="" || !files.ref_seq || !files.soap_result ) {
		// These are compulsory parameters
		usage();
	}
	if ( ! para->glf_format ) {
		// Normal SOAPsnp tab-delimited text format
		files.consensus.clear();
		files.consensus.open(consensus_name.c_str());
		if( ! files.consensus ) {
			cerr<<"Cannot creat file:" <<consensus_name<<endl;
			exit(1);
		}
		files.consensus.clear();
	}
	else {
		// SOAPsnp-defined GLF and baseinfo format
		files.consensus.clear();
		files.consensus.open(consensus_name.c_str(), ios::binary);
		if(!files.consensus) {
			cerr<<"Cannot creat file:" <<consensus_name<<endl;
			exit(1);
		}
		files.consensus.clear();

		files.baseinfo.clear();
		string baseinfo_name = consensus_name + ".baseinfo";
		files.baseinfo.open(baseinfo_name.c_str());
		if(!files.baseinfo) {
			cerr<<"Cannot creat file:" <<baseinfo_name<<endl;
			exit(1);
		}
		files.baseinfo.clear();

		files.o_region.clear();
		string o_region_name = consensus_name + ".index";
		files.o_region.open(o_region_name.c_str());
		if(!files.o_region) {
			cerr<<"Cannot creat file:" <<o_region_name<<endl;
			exit(1);
		}
		files.o_region.clear();
	}

	//Read the chromosomes into memory
	Genome * genome = new Genome(files.ref_seq, files.dbsnp);
	files.ref_seq.close();
	files.dbsnp.close();
	clog<<"Reading Chromosome and dbSNP information Done."<<endl;
	if(para->region_only && files.region) {
		genome->read_region(files.region, para);
		clog<<"Read target region done."<<endl;
	}

	Prob_matrix * mat = new Prob_matrix;
	if( ! is_matrix_in) {
		//Read the soap result and give the calibration matrix
		mat->matrix_gen(files.soap_result, para, genome);
		if (files.matrix_file) {
			mat->matrix_write(files.matrix_file, para);
		}
	}
	else {
		mat->matrix_read(files.matrix_file, para);
	}
	files.matrix_file.close();
	clog<<"Correction Matrix Done!"<<endl;
	mat->prior_gen(para);
	mat->rank_table_gen();
	Call_win info(para->read_length);
	info.initialize(0);
	//Call the consensus
	files.soap_result.close();
	files.soap_result.clear();
	files.soap_result.open(alignment_name.c_str());
	files.soap_result.clear();
	info.soap2cns(files.soap_result, files.consensus, files.baseinfo, genome, mat, para);
	files.soap_result.close();
	files.consensus.close();
	cerr<<"Consensus Done!"<<endl;
	return 0;
}

