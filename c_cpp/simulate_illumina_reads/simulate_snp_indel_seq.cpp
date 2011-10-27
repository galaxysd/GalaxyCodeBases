#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <map>
#include "simulate.h"
#include "simulate_snp_indel_seq.h"
#include "gzstream.h"

using namespace std;

string input;
float hetersnp_rate=0;
float heterindel_rate=0;
int output_type = 1;
ofstream outfile;
ogzstream gz_outfile;
string output = "ref_sequence";

const char *MAKE_TIME="2011-04-22";
const char *VERSION="2.0";
const char *AUTHOR="yuanjianying";
const char *CONTACT="yuanjianying@genomics.org.cn";

void Usage(){
	cout<<"Description:"<<endl;
	cout<<endl<<"\tIt is a program for simulating heterozygosis SNP and heterozygosis Indel in diploid,"<<endl;
	cout<<"\tinput sequence must be set ,because there is not default value."<<endl;
	cout<<endl<<"Program:simulate_snp_indel_seq"<<endl;
	cout<<"\tCompile Data: "<<MAKE_TIME<<endl;
	cout<<"\tAuthor: "<<AUTHOR<<endl;
	cout<<"\tVersion: "<<VERSION<<endl;
	cout<<"\tContact: "<<CONTACT<<endl;
	cout<<endl<<"Usage:\tsimulate_snp_indel_seq [options]"<<endl;
	cout<<"\t-i	<string>  input,input reference genome sequence *.fa/*.fa.gz"<<endl;
	cout<<"\t-s	<float>   heterSNP_rate,set the heterozygous SNP rate of the diploid genome,default:0"<<endl;
	cout<<"\t-d	<float>   heterIndel_rate,set the heterozygous indel rate of the diploid genome,default:0"<<endl;
	cout<<"\t-c	<int>     set ouput file type, 0:text, 1:*.gz, default:1"<<endl;
	cout<<"\t-o	<string>  output,output file prefix default: ref_sequence"<<endl;
	cout<<"\t-h	          help,output help infomation"<<endl;
	cout<<endl<<"Example:"<<endl;
	cout<<endl<<"\t1. ./simulate_snp_indel_seq -i ref_sequence.fa -s 0.001 -d 0.0001 -o ref_sequence >simulate_snp_indel.out 2>simulate_snp_indel.err"<<endl;
	exit(-1);
}

void Getopt(int argc,char *argv[]){
	int c;
	while ((c=getopt(argc,argv,"i:s:d:o:c:h"))!=-1)
	{
		switch(c){
			case 'i': input=optarg;break;
			case 's': hetersnp_rate=atof(optarg);break;
			case 'd': heterindel_rate=atof(optarg);break;
			case 'c': output_type=atoi(optarg);break;
			case 'o': output=optarg;break;
			case 'h': Usage();break;
			default: Usage();
		}
	}
}

int main(int argc, char *argv[])
{
	time_t time_start, time_end;
	time_start = time(NULL);
	
	if (argc==1)
	{
		Usage();
	}
	Getopt(argc,argv);
	
	if(hetersnp_rate < 0 || hetersnp_rate > 1){cerr<<"Error: hetersnp_rate should be set in 0~1 "<<endl;exit(-1);}
	if(heterindel_rate < 0 || heterindel_rate > 1){cerr<<"Error: heterindel_rate should be set in 0~1 "<<endl;exit(-1);}
	if(hetersnp_rate == 0 && heterindel_rate == 0){cerr<<"Please input hetersnp_rate or heterindel_rate!"<<endl; exit(-1);}
		
	igzstream infile;
	infile.open(input.c_str());
	if (!infile)
	{
		cerr<<"Error:unable to open input file:"<<input<<endl;
		exit(-1);
	}
	string snp_output = output+"_snp.lis";
	string indel_output = output+"_indel.lis";
	if(hetersnp_rate != 0){output = output+".snp";}
	if(heterindel_rate != 0){output = output+".indel";}
	if(!output_type){
		output = output + ".fa";
	}else{
		output = output + ".fa.gz";
	}
	
//	ofstream outfile;
	if(!output_type){
		outfile.open(output.c_str());
  	if (!outfile)
  	{
  		cerr<<"Error:unable to open output file:"<<output<<endl;
  		exit(-1);
  	}
	}else{
		gz_outfile.open(output.c_str());
  	if (!gz_outfile)
  	{
  		cerr<<"Error:unable to open output file:"<<output<<endl;
  		exit(-1);
  	}
	}

	ofstream snp,indel;

	if(hetersnp_rate>0){
  	snp.open(snp_output.c_str());
  	if(!snp)
  	{
  		cerr<<"Error:unable to open output file:"<<snp_output<<endl;
  		exit(-1);
  	}
	}
	
	if(heterindel_rate>0){
  	indel.open(indel_output.c_str());
  	if(!indel)
  	{
  		cerr<<"Error:unable to open output file:"<<indel_output<<endl;
  		exit(-1);
  	}
	}
	
//	Get_raw_genome(infile,outfile,snp,indel);
	Get_raw_genome(infile,snp,indel);
	
	time_end = time(NULL);
	cerr<<"All done! Run time: "<<time_end-time_start<<"s."<<endl;
	
	return 0;
}

//get raw genome sequence.
void Get_raw_genome(igzstream &inf, ofstream &snp, ofstream &indel)
{
	string line,id,id_line,seq;
	while (getline(inf,line,'\n'))
	{
		if (line[0]=='>')
		{
			if (seq!="")
			{	
				cerr<<"Have finished reading scaffold "<<id<<endl;
//				simulate_snp_indel_seq(id_line,id,seq,outfile,snp,indel);
				simulate_snp_indel_seq(id_line,id,seq,snp,indel);
				seq="";
			}
			id_line = line;
			line.erase(0,1);
//			id=line;
			int pos=line.find(" ");
			line=line.substr(0,pos);
			id=line;
		}else{
			seq+=line;
		}		
	}
	cerr<<"Have finished reading scaffold "<<id<<endl;
//	simulate_snp_indel_seq(id_line,id,seq,outfile,snp,indel);
	simulate_snp_indel_seq(id_line,id,seq,snp,indel);
}

//add snp and indel in raw seqence, and output result sequence.
void simulate_snp_indel_seq(string id_line,string id,string &sequ, ofstream &snp,ofstream &indel)
{
	srand((unsigned)time(NULL));
	string sequence;
	//convert lower case to upper case 
	for (int i=0;i<sequ.size();i++)
	{
		sequence.push_back(Bases[alphabet[sequ[i]]]);
	}

	if (hetersnp_rate>0 || heterindel_rate>0) //heterozygous SNP and heterozygous indel exists in diploid
	{
		if (hetersnp_rate>0)
		{
			cerr<<"Begin to simulate snp"<<endl;
			sequence=Get_snp(sequence,snp,id,hetersnp_rate);
			cerr<<"Finish to simulate snp"<<endl;
		}
		if (heterindel_rate>0)
		{
			cerr<<"Begin to simulate indel"<<endl;
			sequence=Get_indel(sequence,indel,id,heterindel_rate);
			cerr<<"Finish to simulate indel"<<endl;
		}
		if(!output_type){
			outfile<<id_line<<endl
  			<<sequence<<endl;
		}else{
			gz_outfile<<id_line<<endl
  			<<sequence<<endl;
		}
//		outfile<<id_line<<endl
//			<<sequence<<endl;
		cerr<<"Finish simulate "<<id_line<<endl;
	}
}


