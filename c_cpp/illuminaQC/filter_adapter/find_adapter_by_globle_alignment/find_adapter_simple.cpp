//author: Fan Wei, email: fanw@genomics.org.cn, date: 2007-9-22
//function: find adapter with only one specified insert size
//find_adapter_simple是find_adapter的简单版

#include "common.h"
#include "locate_adapter.h"
#include "gzstream.h"

using namespace std;

string adapter_type = "gDNA-3"; // default adapter type 
int insert_size = 0;
int verbose = 0;
string adpt1="AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG";//pDNA-3+
string adpt2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";//pDNA-5-

void usage ()
{	
	cout << "\nUsage: find_adapter [options] <infile.fq|gz>\n"
			<< "  -a <str>   set adapter type, gDNA-3 for read1, and gDNA-5 for read2, default=" << adapter_type << "\n"
			<< "  -i <int>   set inserted fragment size, default="<< insert_size << "\n"
			<< "  -v   output verbose progressing information\n"
			<< "  -h   outpput help information\n" 
			<< endl ;
	exit(0);
}



int main ( int argc, char *argv[] )
{	
	//deal with command line options
	int c; //must be int type
	while((c=getopt(argc,argv,"a:i:vh")) != -1) { // :表示有参数值
		switch ( c )
		{	// return optarg, type is char*, use atoi or atof to convert
			case 'a' : adapter_type = optarg; break;
			case 'i' : insert_size = atoi(optarg); break;
			case 'v' : verbose = 1; break;
			case 'h' : usage(); break;
			default  : usage(); 
		}
	}
	if (argc < 2) usage();
	
	string adapter;
	if (adapter_type == "gDNA-3")
	{	adapter = adpt1;
	}
	if (adapter_type == "gDNA-5")
	{	adapter = adpt2;
	}
	
	cerr << "adapter type: " << adapter_type << endl;
	cerr << "adapter seq:  " << adapter << endl;
	cerr << "insert size:  " << insert_size << endl;

	
	//get input file and set output file names
	string seq_file = argv[optind++];  //argv[optind++]顺序指向非option的参数
	string res_file = seq_file + ".adapter." + adapter_type + "." + int2str(insert_size) + ".list.gz";
	string sta_file = seq_file + ".adapter." + adapter_type + "." + int2str(insert_size) + ".stat";
	
	//cout << seq_file << " " << res_file << " " << sta_file << endl;


	igzstream infile ( seq_file.c_str() );
	if ( ! infile )
	{	cerr << "fail to open input file " << seq_file << endl;
	}

	ogzstream resfile ( res_file.c_str() );
	if ( ! resfile )
	{	cerr << "fail to creat result file " << res_file << endl;
	}
	
	ofstream stafile ( sta_file.c_str() );
	if ( ! stafile )
	{	cerr << "fail to creat stat file " << sta_file << endl;
	}

	int total_num = 0; //number of total reads 
	int adapter_num = 0; //number of reads with adapter
	float adapter_rate = 0.0; //rate of reads adapter
	
	resfile << "id    seq    qual    align_len    mis_match    insert_size\n";

	//parse one sequence at a time
	string textLine;
	while ( getline( infile, textLine, '\n' ) )
	{	if ( textLine[0] == '@' )
		{	
			string id, seq, qual;
			vector <string> vecLine;

			split(textLine,vecLine,"@");
			id = vecLine[0];
			
			getline( infile, textLine, '\n' );
			seq = textLine;

			getline( infile, textLine, '\n' );
			getline( infile, textLine, '\n' );
			qual = textLine;
			
			if (seq.size() < 5 || adapter.size() < 5) continue; 

			int align_len=0, mis_match=0;

			if ( locate_adapter(seq,adapter,insert_size,align_len, mis_match) )
			{
				resfile << id << "  " << seq << "  " << qual << "  " << align_len << "  " << mis_match << "  " << insert_size << "\n";
				adapter_num++;
			}

			total_num++;
		} 
	}
	
	adapter_rate = float(adapter_num) / total_num;

	stafile << "total reads    :  " <<  total_num << endl;
	stafile << "adapter reads  :  " << adapter_num << endl;
	stafile << "adapter rate   :  " << float2str(adapter_rate)  << endl;

}


