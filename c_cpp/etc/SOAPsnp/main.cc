/**
  *  SOAPsnp (Short Oligonucleotide Analysis Package for Single Nucleotide Polymorphism)
  *
  *  main.cc
  *
  *  Copyright (C) 2008, BGI Shenzhen.
  *
  *
  */

#include "soap_snp.h"
#include <getopt.h>

using namespace std;

const char *PROGRAM = "SOAPsnp (Short Oligonucleotide Analysis Package for Single Nucleotide Polymorphism)";
const char *AUTHOR = "BGI Shenzhen";
const char *VERSION = "1.00";
const char *CONTACT = "soap@genomics.org.cn";

int usage() {
	cerr<<"\nProgram: "<<PROGRAM<<endl;
	cerr<<"Copyright (C) 2008, BGI Shenzhen."<<endl;
	cerr<<"Author:  BGI Shenzhen"<<endl;
	cerr<<"Version: "<<VERSION<<endl;
	cerr<<"Contact: soap@genomics.org.cn\n"<<endl;
	cerr<<"Compulsory Parameters:"<<endl;
	cerr<<"-i <FILE> Input SORTED soap Result"<<endl;
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
	cerr<<"-n Enable binomial probability calculation to give HET for better accuracy. [Off]"<<endl;
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
	cerr<<"-h Display this help\n"<<endl;
	exit(1);
	return 0;
}

int readme() {
	return usage();
}

int main ( int argc, char * argv[]) {
	// This part is the default values of all parameters
	Parameter * para = new Parameter;
	std::string alignment_name, consensus_name;
	bool is_matrix_in = false; // Generate the matrix or just read it?
	int c;
	Files files;
	while((c=getopt(argc,argv,"i:d:o:z:g:p:r:e:ts:2a:b:j:k:unmqM:I:L:Q:S:F:E:T:h")) != -1) {
		switch(c) {
			case 'i':
			{
				// Soap Alignment Result
				files.soap_result.open(optarg);
				files.soap_result.clear();
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
				files.ref_seq.open(optarg);
				files.ref_seq.clear();
				if( ! files.ref_seq) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					exit(1);
				}

				break;
			}
			case 'o':
			{
				files.consensus.open(optarg);
				files.consensus.clear();
				if( ! files.consensus ) {
					cerr<<"Cannot creat file:" <<optarg <<endl;
					exit(1);
				}
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
				files.dbsnp.open(optarg);
				files.dbsnp.clear();
				if( ! files.ref_seq) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					exit(1);
				}
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
				if (files.matrix_file) {
					cerr<<"Plz do not use -M and -I simutaneously"<<endl;
					exit(1);
				}
				// Output the calibration matrix
				files.matrix_file.open(optarg, fstream::out);
				files.matrix_file.clear();
				if( ! files.matrix_file) {
					cerr<<"Cannot creat file :"<<optarg<<endl;
					exit(1);
				}
				break;
			}
			case 'I':
			{
				if (files.matrix_file) {
					cerr<<"Plz do not use -M and -I simutaneously."<<endl;
					exit(1);
				}
				// Input the calibration matrix
				files.matrix_file.open(optarg, fstream::in);
				files.matrix_file.clear();
				if( ! files.matrix_file) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					exit(1);
				}
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
	if( !files.consensus || !files.ref_seq || !files.soap_result ) {
		// These are compulsory parameters
		usage();
	}
	//Read the chromosomes into memory
	Genome * genome = new Genome(files.ref_seq, files.dbsnp);
	files.ref_seq.close();
	files.dbsnp.close();
	clog<<"Reading Chromosome and dbSNP information Done."<<endl;
	if(para->region_only && files.region) {
		genome->read_region(files.region);
		clog<<"Read target region done."<<endl;
	}
	if(para->glf_format) { // GLF or GPF
		files.consensus.close();
		files.consensus.clear();
		files.consensus.open(consensus_name.c_str(), ios::binary);
		if(!files.consensus) {
			cerr<<"Cannot write result ot the specified output file."<<endl;
			exit(255);
		}
		if (1==para->glf_format) {
			files.consensus<<'g'<<'l'<<'f';
		}
		else if (2==para->glf_format) {
			files.consensus<<'g'<<'p'<<'f';
		}
		int major_ver = 0;
		int minor_ver = 0;
		files.consensus.write(reinterpret_cast<char*>(&major_ver), sizeof(major_ver));
		files.consensus.write(reinterpret_cast<char*>(&minor_ver), sizeof(minor_ver));
		if(!files.consensus.good()) {
			cerr<<"Broken ofstream after version."<<endl;
			exit(255);
		}
		std::string temp("");
		for(std::string::iterator iter=para->glf_header.begin();iter!=para->glf_header.end(); iter++) {
			if (':'==(*iter)) {
				int type_len(temp.size()+1);
				files.consensus.write(reinterpret_cast<char*>(&type_len), sizeof(type_len));
				files.consensus.write(temp.c_str(), temp.size()+1)<<flush;
				temp = "";
			}
			else {
				temp+=(*iter);
			}
		}
		if(!files.consensus.good()) {
			cerr<<"Broken ofstream after tags."<<endl;
			exit(255);
		}
		if(temp != "") {
			int type_len(temp.size()+1);
			files.consensus.write(reinterpret_cast<char*>(&type_len), sizeof(type_len));
			files.consensus.write(temp.c_str(), temp.size()+1)<<flush;
			temp = "";
		}
		int temp_int(12);
		files.consensus.write(reinterpret_cast<char*>(&temp_int), sizeof(temp_int));
		files.consensus.write("CHROMOSOMES", 12);
		temp_int = genome->chromosomes.size();
		files.consensus.write(reinterpret_cast<char*>(&temp_int), sizeof(temp_int));
		files.consensus<<flush;
		if(!files.consensus.good()) {
			cerr<<"Broken ofstream after writting header."<<endl;
			exit(255);
		}
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
	info.soap2cns(files.soap_result, files.consensus, genome, mat, para);
	files.soap_result.close();
	files.consensus.close();
	cerr<<"Consensus Done!"<<endl;
	return 0;
}

