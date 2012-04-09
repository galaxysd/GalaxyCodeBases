//author: Fan Wei, email: fanw@genomics.org.cn, date: 2007-9-22
//function: find adapter with multiple insert size
//本程序相当于只做正链的序列比对，但没有用动态规划算法，而是挨个碱基挨个位置的比较。

#include "common.h"
#include "locate_adapter.h"
#include "gzstream.h"

using namespace std;

string adapter_type = "gDNA-3"; // default adapter type 
int insert_size = 0;
int verbose = 0;
string adpt1="AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG";//pDNA-3+
string adpt2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";//pDNA-5-


void usage ();
void gnuplot_drawing (map <int,int> &adapter_num, string file_core, int align_start, int align_end, int total_num);

int main ( int argc, char *argv[] )
{
	//deal with command line options
	int c; //must be int type
	while((c=getopt(argc,argv,"a:i:vh")) != -1) { // :表示有参数值
		switch ( c )
		{	// return optarg, type is char*, use atoi or atof to convert
			case 'a' : adapter_type = optarg; break;
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
	

	//get input file and set output file names
	string seq_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string res_file = seq_file + ".adapter." + adapter_type + "." + ".list.gz";
	string sta_file = seq_file + ".adapter." + adapter_type + "." + ".stat";
	

	int total_num = 0; //number of total reads
	int adapter_all = 0; //number of reads with adapter of all insert sizes
	float adapter_rate = 0.0; //rate of reads with adapter
	map <int,int> adapter_num; //insert size, adapter reads
	int max_seq_len = 0; //the longest read in a fastq file
	//cout << seq_file << " " << res_file << " " << sta_file << endl;

	igzstream infile ( seq_file.c_str() );
	if ( ! infile )
	{	cerr << "fail to open input file" << seq_file << endl;
	}

	ogzstream resfile ( res_file.c_str() );
	if ( ! resfile )
	{	cerr << "fail to creat result file" << res_file << endl;
	}
	
	ofstream stafile ( sta_file.c_str() );
	if ( ! stafile )
	{	cerr << "fail to creat stat file" << sta_file << endl;
	}
	
	resfile << "id    seq    qual    align_len    mis_match　 insert_size\n";

	//parse one sequence at a time
	string id,seq,qid,qual;
	while ( getline( infile, id, '\n' ) )
	{	if ( id[0] == '@' )
		{				
			getline( infile, seq, '\n' );
			getline( infile, qid, '\n' );
			getline( infile, qual, '\n' );

			int align_len=0, mis_match=0;
			int is_adapter = 0;
			
			if (seq.size() < 5 || adapter.size() < 5) continue; 
			if (max_seq_len < seq.size()) max_seq_len = seq.size();

			//first search insert size from 0 to positive number
			for (int insert_size=0; insert_size<seq.size()-5; insert_size++)
			{	is_adapter = locate_adapter(seq,adapter,insert_size,align_len, mis_match);
				if ( is_adapter )
				{	
					if (adapter_num.count(insert_size))
					{	
						adapter_num[insert_size]++;
					}else{
						adapter_num[insert_size] = 1;
					}
					adapter_all++;
					
					resfile << id << "  " << seq << "  " << qual << "  " << align_len << "  " << mis_match << "  " << insert_size << endl;
					break;
				}
			}


			//如果在0及正数insert size里面没找到adapter
			//then search insert size in negative number
			if (! is_adapter)
			{	for (int insert_size=-1; abs(insert_size) < adapter.size()-5; insert_size--)
				{	
					is_adapter = locate_adapter(seq,adapter,insert_size,align_len, mis_match);
					if ( is_adapter )
					{	
						if (adapter_num.count(insert_size))
						{	
							adapter_num[insert_size]++;
						}else{
							adapter_num[insert_size] = 1;
						}
						adapter_all++;

						resfile << id << "  " << seq << "  " << qual << "  " << align_len << "  " << mis_match << "  " << insert_size << endl;

						break;
					}
				}
			}

			total_num++;
		} 
	}

	//static the results
	adapter_rate = float(adapter_all) / total_num;
	int empty_num = (adapter_num.count(0)) ? adapter_num[0] : 0;
	float empty_rate = float(empty_num) / total_num;

	stafile << "total reads   :  " << total_num << "\n";
	stafile << "adapter reads :  " << adapter_all << "\n";
	stafile << "adapter rate  :  " << float2str(adapter_rate) << "\n";
	stafile << "empty reads   :  " << empty_num << "\n";
	stafile << "empty rate    :  " << float2str(empty_rate) << "\n\n";
	
	stafile << "insert_size\tadapter_reads\tadapter_rate\n"; 
	
	int align_start = 5 - adapter.size();
	int align_end = max_seq_len - 5;
	for (int i=align_start; i<=align_end; i++)
	{	
		int insert_size_num = (adapter_num.count(i)) ? adapter_num[i] : 0;
		float insert_size_rate = float(insert_size_num) / total_num;
		stafile << i << "\t" << insert_size_num << "\t" << float2str(insert_size_rate) << "\n";
	}
	
	//draw a distribution figure for adapter rate
	//gnuplot_drawing(adapter_num, sta_file, align_start, align_end, total_num);

}


void gnuplot_drawing (map <int,int> &adapter_num, string file_core, int align_start, int align_end, int total_num)
{	
	string fig_file = file_core + ".png";
	string dat_file = file_core + ".dat";
	string gnuplot_file = file_core + ".gnuplot";
	
	ofstream datfile ( dat_file.c_str() );
	if ( ! datfile )
	{	cerr << "fail to creat result file" << dat_file << endl;
	}
	ofstream gnufile ( gnuplot_file.c_str() );
	if ( ! gnufile )
	{	cerr << "fail to creat result file" << gnuplot_file << endl;
	}
	
	//generate datfile for gnuplot
	for (int i=align_start; i<=align_end; i++)
	{	
		int insert_size_num = (adapter_num.count(i)) ? adapter_num[i] : 0;
		float insert_size_rate = float(insert_size_num) / total_num;
		datfile << i << "\t" << float2str(insert_size_rate) << endl;
	}
	
	//generate gnuplot shell file
	gnufile << "set title \"Adapter Pollution\"\n"
		 << "set xlabel \"insert size\"\n"
		 << "set ylabel \"adapter rate\"\n"
		 << "set terminal png\n"         
		 << "set output \"" << fig_file << "\"\n"
		 << "plot \"" << dat_file 
		 << "\" using 1:2 title \"\" with linespoints;\n"
		 << endl;

	string command_line = "/usr/local/bin/gnuplot " + gnuplot_file + "; ";
	command_line += "rm " + gnuplot_file  + "; ";
	command_line += "rm " + dat_file  + "; ";
	system(command_line.c_str());

}


void usage ()
{	
	cout << "\nUsage: find_adapter [options] <infile.fq|gz> \n"
			<< "  -a <str>   set adapter type, gDNA-3 for read1, and gDNA-5 for read2, default=" << adapter_type << "\n"
			<< "  -v         output verbose progressing information\n"
			<< "  -h         outpput help information\n" 
			<< endl ;
	exit(1);
}
