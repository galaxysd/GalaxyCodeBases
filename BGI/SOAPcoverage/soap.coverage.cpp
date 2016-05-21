/*
	Enchantez!

	Soap.Coverage

	Author: Aqua Lo
	E-mail: luoruibang@genomics.org.cn

	This utility if a component of SOAP.
*/

#define _REENTRANT
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <cassert>
#include <limits>
#include <pthread.h>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>
#include <boost/regex.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include "gzstream.h"
#include "threadmanager.h"


//Namespace Declaration
using namespace std;
using namespace boost;

//Constant
#define DEPTH_SINGLE_LENGTH_PER_ROW 100		//Columns per row in -depthsingle
#define WINDOW_SIZE_LOWER_LIMIT 100			//The lowest threshold of window_size
//Debug Preparations
#define DEBUG false;
#define DB(...) if(DEBUG) cerr<<"Debug:"<<'\t'<<__FILE__<<'\t'<<__LINE__<<'\t'<<(__VA_ARGS__)<<end;


//Input Method
typedef igzstream imethod;

//Program Label
#define PROGRAM_NAME "SOAP.coverage"
#define PROGRAM_COMPILE_DATE __DATE__
#define PROGRAM_COMPILE_TIME __TIME__

// Define version number
#define VERSION_MAJOR 2
#define VERSION_MINOR 7
#define VERSION_PATCHLEVEL 9
#define VERSION_STRING "2.7.9"

// Macros dealing with VERSION
#define VERSION_NUM(a,b,c) (((a) << 16L) | ((b) << 8) | (c))
#define VERSION \
		VERSION_NUM(VERSION_MAJOR, VERSION_MINOR, VERSION_PATCHLEVEL)

#define README \
	"This utility can calculate sequencing coverage or physical coverage as well as duplication rate\n"\
	"and details of specific block for each segments and whole genome by using SOAP, Blat, Blast, BlastZ,\n"\
	"mummer and MAQ aligement results with multi-thread. Gzip file supported.\n"

#define USAGE \
	"Parameters:\n"\
	"  -cvg or -phy or -tag    Selector for sequencing coverage mode, physical coverage mode or reads tag mode\n"\
	"                          At least and only one should be selected!\n"\
	"  -refsingle [filename]   Input reference fasta file used in SOAP\n"\
	"  -i [soap-file1 soap-file2 ...]\n"\
	"                          Input several soap or soap gziped results by filenames.\n"\
	"  -il [soap-list]         Input several soap or soap gziped results (absolute path!) with a soap-list file\n"\
	"		Caution: Only PE aligned results can be used in physical coverage!\n"\
	"  -il_single [SE aligned results list]\n"\
	"  -il_soap [PE aligned results list]\n"\
	"  -o [file-name]          Results output with details\n"\
	"  -depth [directory]      Output coverage of each bp in seperate files, a directory should be given\n"\
	"  -depthsingle [filename] Output coverage of each bp in a single file (text, fasta like)\n"\
	"  -nospace                No space between plain depthsingle site depths\n"\
	"  -depthsinglebin [fn]    Output coverage of each bp in a single file (Binary mode)\n"\
	"  -addn [filename]        Input N block data for exclusion (marked as 65535 in depthsingle output)\n"\
	"		Input format: <segment_name> <start (leftmost as 1)> <end (half close)>\n"\
	"  -depthinput [filename]  Input previous coverage data (Both Text or Binary) for faster accumulation\n"\
	"  -depthinputsam [filename]\n"\
	"                          Input previous coverage data from Samtools\n"\
	"  -cdsinput [filename]    Input specific block range for calculating coverage\n"\
	"		Input format: <segment_name> <start (leftmost as 1)> <end (half close)>\n"\
	"  -plot [filename] [x-axis lower] [x-axis upper]\n"\
	"                          Output overall distribution of coverage of all segments\n"\
	"  -cdsplot [filename] [x-axis lower] [x-axis upper]\n"\
	"                          Output distribution of coverage of specific blocks\n"\
	"  -cdsdetail [filename]\n"\
	"                          Output coverage of each bp of each specific blocks in a single file\n"\
	"  -window [filename] [length]\n"\
	"                          Output coverage averaged in a [length] long window to [filename]\n"\
	"  -p [num]                Number of processors [Default:4]\n"\
	"  -trim [num]             Exclude [num] bp(s) from head & tail of each segments\n"\
	"  -precisionOffset [num]  Set the precision offset [Default:2]\n"\
	"  -minimumDepth [num]     Set the minimum depth to take into consideration [Default:0]\n"\
	"  -maximumDepth [num]     Set the maximum depth to take into consideration [Default:65535]\n"\
	"\n"\
	"Input format seletors:\n"\
	"  -general [id column] [start position column] [end position column]\n"\
	"  -generallen [id column] [start position column] [len]\n"\
	"                          Specify column's (leftmost 1 based, right half close) for necessary informations\n"\
	"                              Warning: Inputs with /^#/ are recognized as comments and ignored\n"\
	"  -plain                  Input is a three column list\n"\
	"  -sam                    Input is a standard SAM input file\n"\
	"  -pslquery               Input is Blat for alculating query coverage.\n"\
	"  -pslsub                 Input is Blat for calculating subject coverage.\n"\
	"  -maq                    Input is MAQ output file.\n"\
	"  -m8subject              Input is Blast m8 file for calculating subject coverage (reference should be subject).\n"\
	"  -m8query                Input is Blast m8 file for calculating query coverage (reference should be query).\n"\
	"  -mummerquery [limit]    Input mummer result file for calculating query coverage.\n"\
	"  -axtoitg                Input Blastz axt file for calculating target coverage.\n"\
	"  -axtoiq                 Input Blastz axt file for calculating query coverage.\n"\
	"  -qmap                   Input qmap of bisulfite sequencing.\n"\
	"\n"\
	"Special functions:\n"\
	"  -sp [filename_in] [filename_out]\n"\
	"                          Output S/P ratio data for post processing.\n"\
	"       Column:\n"\
	"           ref    start    end    name\n"\
	"  -pesupport [filename_in] [filename_out]\n"\
	"                          Output pair-end reads on specific areas.\n"\
	"       Column:\n"\
	"           ref    start    end    name\n"\
	"  -onlyuniq               Use reads those are uniquely mapped (column 4 in soap == 1).\n"\
	"  -precise                Omit mismatched bp in soap results.\n"\
	"  -nowarning              Cancel all possible warning.\n"\
	"  -nocalc                 Do not perform depth calculation.\n"\
	"  -onlycover              Only output 0 or 1 for coverage calculation.\n"\
	"\n"\
	"Physical Coverage Specified Parameters:\n"\
	"  -duplicate [num]        Exclude duplications, and gives the percentage of duplication. [num]=readlength\n"\
	"  -insertupper [num]      Insert larger than num will be abandon [Default: 15000]\n"\
	"  -insertlower [num]      Insert shorter thab num will be abandon [Default: 0]\n"\
	"\n"\
	"Example:\n"\
	"	1. Calculate several files of SOAP results.\n"\
	"	   soap.coverage -cvg -i *.soap *.single -refsingle human.fa -o result.txt \n"\
	"\n"\
	"	2. Calculate a list of SOAP results, exclude Ns blocks, output depth to\n"\
	"	   a file and plot coverage form depth 0 to 1000.\n"\
	"	   soap.coverage -cvg -il soap.list -refsingle human.fa -o result.txt -depthsingle all.coverage -addn n.gap -plot distribution.txt 0 1000\n"\
	"\n"\
	"	3. Calculate a list of SOAP results, use only uniquely mapped reads, exclude Ns blocks\n"\
	"	   , output depth to a file and plot coverage form depth 0 to 1000.\n"\
	"	   soap.coverage -cvg -il soap.list -onlyuniq -refsingle human.fa -o result.txt -depthsingle all.coverage -addn n.gap -plot distribution.txt 0 1000\n"\
	"\n"\
	"	4. Add new SOAP results to depth(-depthsingle) already calculated &\n"\
	"	   plot all data and specific blocks from depth 0 to 150, with 6 processors.\n"\
	"	   soap.coverage -cvg -depthinput all.coverage -refsingle human.fa -il soap.list -p 6 -o result.txt -cdsinput cds.list -plot distribution.txt 0 150 -cdsplot distribution_cds.txt 0 150\n"\
	"\n"\
	"	5. Calculate physical coverage and duplication rate(read length=44) with\n"\
	"	   insert between (avg-3SD, avg+SD)[avg=197, SD=9], with 8 processors\n"\
	"	   soap.coverage -phy -il soap_without_single.list -refsingle human.fa -p 8 -o result.txt -duplicate 44 -insertlower 170 -insertupper 224\n"\

//Typedef
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unordered_set<string> strsets;
typedef unordered_map<string,int> strmaps;
typedef unordered_map<string,pthread_mutex_t*> semaphore_maps;
typedef vector<ushort> vs;
typedef vector<uint> vi;
typedef vector<pair<string, uint> > vpsu;

//Functions Declaration
void Intro();								//Program Head
void ParaDeal(int&, char**);				//Dealing with parameters
void MakeMaps();							//Segments Shadowing form sets to maps
void BuildMemory();							//Alloc Memory
inline ulong AddTotalBP();					//Use by BuildMemory(), cumulate quantity info in ref.
inline ulong AddTotalBP_refsingle();		//Use by BuildMemory(), cumulate quantity info in refsingle.
void AddCover();							//-cvg
void Physical_AddCover();					//-phy
void CalcCoverage();
void OutputCoverage();						//-depthsingle
void OutputCoverageBin();					//-depthsinglebin
void OutputDistribution();					//-plot
void OutputWindow();						//-window
void AddN();								//-addn
void CalcCDS();								//-cdsplot
void CalcPoint();							//-point
void SP();									//-sp
void SPOutput();
void PESupportOutput();						//-pesupport
void AddCvg();								//-depthinput
void AddCvgFromSamtools();					//-depthinputsam
void AddRefSingle();						//-refsingle
void AddFastat();							//-reffastat
void Output_list();							//Sub of result illustration
void ComfirmAllDigit(const char*);			//Sub of ParaDeal()
inline uint SeperateCigar(string &, vector<uint> &);			//Sub of _AddCover_Sub

//Global vars
bool cvg(false), phy(false), tag(false);
bool depth(false), depthsingle(false), nospace(false), depthsinglebin(false), ref(false), refsingle(false), reffastat(false),depthinput(false),depthinputsam(false),\
		cdsinput(false), cdsplot(false), addon(false), addN(false), plot(false),cdsdetail(false), sam(false),\
		window(false), pslsub(false), pslquery(false), maq(false), m8subject(false), m8query(false), axtoitg(false), axtoiq(false),\
		nowarning(false), sp(false), il_single(false), il_soap(false), pesupport(false), onlyuniq(false),\
		precise(false), gz(false), mummerquery(false), nocalc(false), onlycover(false), plain(false), qmap(false), general(false),\
		generalLen(false);
map<string, ulong> refsingle_marker;
int trim(0), sprd_length(0), mummer_limit(0);
uint length(0);
vector<string> iPos, isPos, ipPos;
string oPos, rPos, dPos, cPos, cdsPos, cdsPos_out, cdsdetailPos, pointPos, pointOut, nPos, plotPos, windowPos_out, spIn, spOut;
uint plot_upper(0), plot_lower(0), cds_plot_upper(0), cds_plot_lower(0);
uint window_size(0);
ulong totalBP(0),countBP(0);
uint insert_Limit(15000);
uint insert_LowerLimit(0);
uint duplicate(0), useful_data(0), useless_data(0);
uint precisionOffset(6);
ushort minimumDepth(0);
ushort maximumDepth(65535);
string temp;								//Quick inter-funtions messenger for BuildMemory().
ofstream OutData;							//Global output handler
strsets chr;
strmaps chrmaps;
semaphore_maps chr_sema;
vector<vs> vector_vs;
vector<vi> sp_vvs;
vpsu sp_sui_single;
vpsu sp_sui_soap;
vector<int> distribution(65536,0);
uint intmax(numeric_limits<int>::max());
uint idCol(1), startCol(2), endCol(3);

//Multi-Thread global var
int NUM_THREADS(4);
void* _AddCover_general(void*);
void* _AddCover_plain(void*);
void* _AddCover_sam(void*);
void* _AddCover_Sub(void*);
void* _AddCover_pslsub(void*);
void* _AddCover_pslquery(void*);
void* _AddCover_maq(void*);
void* _AddCover_m8subject(void*);
void* _AddCover_m8query(void*);
void* _AddCover_mummerquery(void*);
void* _AddCover_axtoitg(void*);
void* _AddCover_axtoiq(void*);
void* _AddCover_qmap(void*);
void* _Physical_AddCover_Sub(void*);

struct param_struct
{
	string pos;
};

int main(int argc, char ** argv)
{
	timer elapsedtime;

	Intro();

	ParaDeal(argc, argv);

	MakeMaps();

	BuildMemory();

	if(sp) SP();

	if(depthinput) AddCvg();
	if(depthinputsam) AddCvgFromSamtools();
	if(addN) AddN();
	if(addon && cvg) AddCover();
	if(addon && phy) Physical_AddCover();

	if(sp && !pesupport) SPOutput();
	if(sp && pesupport) PESupportOutput();

	if(!nocalc) CalcCoverage();

	OutData<<"Overall:"<<endl
			<<"Total:"<<totalBP<<endl
			<<"Covered:"<<countBP<<endl
			<<"Percentage:"<<(double)countBP/totalBP*100<<endl<<endl;

	if(duplicate)
		OutData<<"Duplication rate:"<<(double)useless_data/(useful_data+useless_data)*100<<'%'<<endl;

	if(plot) OutputDistribution();

	if(window) OutputWindow();

	if(depth || depthsingle) OutputCoverage();

	if(depthsinglebin) OutputCoverageBin();

	if(cdsinput) CalcCDS();
	
	if(!pointPos.empty() && !pointOut.empty()) CalcPoint();

	cout<<"Overall:"<<endl
		<<"Total:"<<totalBP<<endl
		<<"Covered:"<<countBP<<endl
		<<"Percentage:"<<(double)countBP/totalBP*100<<endl;

	if(duplicate)
		cout<<"Duplication rate:"<<(double)useless_data/(useful_data+useless_data)*100<<'%'<<endl;

	Output_list();

	cout<<(elapsedtime.elapsed())<<"s Elapsed!"<<endl;

	return EXIT_SUCCESS;
}

void Intro()
{
	cerr<<endl;
	cerr<<PROGRAM_NAME<<endl;
	cerr<<"    Version: "<<VERSION_STRING<<endl;
	cerr<<"Complied at: "<<PROGRAM_COMPILE_DATE<<' '<<PROGRAM_COMPILE_TIME<<endl;
	cerr<<"     Author: RuiBang Luo"<<endl;
	cerr<<"     E-mail: luoruibang@genomics.org.cn"<<endl<<endl;
	cerr<<README<<endl;
}

void OutputWindow()
{
	if(window_size < WINDOW_SIZE_LOWER_LIMIT)
	{
		cerr<<"User defined window length is lower than the lowest threshold "<<WINDOW_SIZE_LOWER_LIMIT<<endl
			<<"Window length will set to "<<WINDOW_SIZE_LOWER_LIMIT<<endl;
		window_size = WINDOW_SIZE_LOWER_LIMIT;
	}
	progress_display ptcount(chrmaps.size(), cout, "Calculating Coverage in window...\n");
	ofstream O(windowPos_out.c_str());
	if(!O)
	{
		cerr<<"Error opening window output file: "<<windowPos_out<<endl;
		return;
	}
	strmaps::iterator mpos;
	for(mpos = chrmaps.begin(); mpos != chrmaps.end(); ++mpos)
	{
		++ptcount;
		for(register int i(1),limit(vector_vs[mpos->second].size()); i<limit;)
		{
			register int BPTotal(0);
			register int BPCount(0);
			register int j( (i+window_size) < limit ? i+window_size : limit);
			if(j - i < window_size * 0.5)
				break;
			for(; i<j ; ++i)
			{
				if(vector_vs[mpos->second][i]==65535)
					continue;
				if(vector_vs[mpos->second][i] < minimumDepth)
					continue;
				if(vector_vs[mpos->second][i] > maximumDepth)
					continue;
				BPTotal += vector_vs[mpos->second][i];
				++BPCount;
			}
			O<<mpos->first<<'\t'<<i-window_size<<'\t'<<j-1<<'\t'<<BPCount<<'\t'<<(double)BPTotal / BPCount<<endl;
		}
	}
}

void Physical_AddCover()
{
	cout<<"Summarizing Physical Coverage..."<<endl;
	CThreadManager threadMgr(NUM_THREADS);
    threadMgr.SetTimer(0.5);
	vector<param_struct> vparam(iPos.size());
	for(register uint i(0); i<iPos.size(); ++i)
	{
		vparam[i].pos = iPos[i];
		threadMgr.AddThread( _Physical_AddCover_Sub, (void*)&vparam[i]);
	}
	threadMgr.Run();
}

void* _Physical_AddCover_Sub(void* param)
{
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos;
	if((pParam->pos).find(".single") != string::npos)
	{
		cout<<" Omitted!"<<endl;
		return (void *) NULL;
	}
	else
		cout<<endl;
	string temp;
	string ctemp, ctemp1, ctemp2, ctemp3, ctemp4, ctemp5, ctemp6, ctemp_pos;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	imethod InData_Sub;
	InData_Sub.open(pParam->pos.c_str());
	//Var for duplication
	multimap<uint,uint> dup_atom;
	vector<multimap<uint,uint> > dup_vector(chrmaps.size());
	//multimap<uint,uint>::iterator p_dup_atom;
	pair<multimap<uint,uint>::iterator , multimap<uint,uint>::iterator > p_dup_atom_pair;
	//***************************************
	while(InData_Sub)
		{
			register uint end(0),start(0),read_length_neg(0),read_length_pos(0), start_dup(0), end_dup(0);
			strmaps::iterator mpos;
			for(;InData_Sub;)
			{
				++useful_data;
				bool continue_mark(false);
				//- 6 for length
				InData_Sub>>ctemp;InData_Sub>>ctemp;InData_Sub>>ctemp;InData_Sub>>ctemp;InData_Sub>>ctemp;InData_Sub>>ctemp5;
				//- 7 for strand
				InData_Sub>>ctemp_pos;
				//- 8 for chr
				InData_Sub>>ctemp1;
				//- 9 for position
				InData_Sub>>ctemp2;InData_Sub.ignore(100000,'\n');
				//+ 6 for length
				InData_Sub>>ctemp;InData_Sub>>ctemp;InData_Sub>>ctemp;InData_Sub>>ctemp;InData_Sub>>ctemp;InData_Sub>>ctemp6;
				//+ 8 for chr
				InData_Sub>>ctemp;InData_Sub>>ctemp3;
				//+ 9 for position
				InData_Sub>>ctemp4;InData_Sub.ignore(100000,'\n');

				//Judgement
				if(ctemp1 != ctemp3)
				{
					getline(InData_Sub, ctemp, '\n');
					cerr<<"Unmatched chromosome context!"<<"\t"<<ctemp1<<"\t"<<ctemp3<<endl;
					continue;
				}
				if(atoi(ctemp2.c_str()) > atoi(ctemp4.c_str()))
				{
					swap(ctemp5, ctemp6);
					swap(ctemp2, ctemp4);
				}
				read_length_neg = atoi(ctemp6.c_str());
				read_length_pos = atoi(ctemp5.c_str());
				if(read_length_neg == 0 || read_length_pos == 0 || ctemp2.size() > 9 || ctemp4.size() > 9)
				{
					cerr<<"One row bypassed!"<<endl;
					continue;
				}
				start = atoi(ctemp2.c_str());
				end = atoi(ctemp4.c_str()) + read_length_neg - 1;
				if(end<=start)
				{
					cerr<<"end<=start"<<endl;
					continue;
				}
				if((end-start)>insert_Limit || (end-start)<insert_LowerLimit) continue;
				temp = ctemp1;
				mpos = chrmaps.find(temp); semapos = chr_sema.find(temp);
				if(mpos == chrmaps.end()) continue;
				if (duplicate)
				{
					start_dup = start; end_dup = end;
					if(read_length_neg < duplicate)
					{
						start_dup = start_dup - (duplicate - read_length_neg -1);
					}
					if(read_length_pos < duplicate)
					{
						end_dup = atoi(ctemp4.c_str()) + duplicate - 2;
					}
					p_dup_atom_pair = dup_vector[mpos->second].equal_range(start_dup);
					if (p_dup_atom_pair.first == p_dup_atom_pair.second)
					{
						dup_vector[mpos->second].insert(make_pair<uint, uint>(start_dup, end_dup));
					}
					else
					{
						for(;p_dup_atom_pair.first != p_dup_atom_pair.second; ++p_dup_atom_pair.first)
						{
							if(p_dup_atom_pair.first->second == end_dup)
							{
								++useless_data;
								continue_mark = true;
								break;
							}
						}
						if (!continue_mark) dup_vector[mpos->second].insert(make_pair<uint, uint>(start_dup, end_dup));
					}
				}
				if (continue_mark) continue;

				break;
			}
			if(start >= vector_vs[mpos->second].size()) continue;
			if(end >= vector_vs[mpos->second].size()) end = vector_vs[mpos->second].size() - 1;
			if(pesupport)
			{
				int cache(0);
				for(register size_t i(start); i <= end; ++i)
				{
					if(sp_vvs[mpos->second][i] == 0)
						continue;
					if(sp_vvs[mpos->second][i] == cache)
						continue;
					cache = sp_vvs[mpos->second][i];
					++(sp_sui_soap[sp_vvs[mpos->second][i]].second);
					sp_sui_soap[sp_vvs[mpos->second][i]].first += ("\t" + lexical_cast<string>(end - start + 1 - read_length_neg - read_length_pos));
				}
			}

				if(!InData_Sub) break;
				if(start >= vector_vs[mpos->second].size()) continue;
				if(end >= vector_vs[mpos->second].size()) end = vector_vs[mpos->second].size()-1;
				for(register uint j(start),k(end); j<=k; ++j)
					++(vector_vs[mpos->second][j]);
		}
}

void OutputDistribution()
{
	ofstream O(plotPos.c_str());
	if(!O)
	{
		cerr<<"Error opening plot output file: "<<plotPos<<endl;
		exit(EXIT_FAILURE);
	}
	for(register int i(plot_lower); i<=plot_upper; ++i)
		//if(distribution[i])
			O<<setw(5)<<i<<setw(12)<<distribution[i]<<endl;
}

void OutputCoverage()
{
	progress_display ptcount(chrmaps.size(), cout, "Output Coverage to files (Text)...\n");
	string temp;
	strmaps::iterator mpos;
	string delim = nospace ? "": " ";
	if (depth)
	{
		for(mpos = chrmaps.begin(); mpos != chrmaps.end(); ++mpos)
		{
			++ptcount;
			temp = dPos + '/'+ mpos->first + ".coverage";
			ofstream depthData(temp.c_str());
			depthData<<'>'<<mpos->first<<"\n";

			//Determine Length
			length = DEPTH_SINGLE_LENGTH_PER_ROW;

			if(onlycover)
			{
				for(register unsigned long i(1),limit(vector_vs[mpos->second].size()); i<limit; ++i)
				{
					if(vector_vs[mpos->second][i])
						depthData<<1<<delim;
					else
						depthData<<0<<delim;
					if (i && i%length == 0) depthData<<'\n';
				}
			}
			else
			{
				for(register unsigned long i(1),limit(vector_vs[mpos->second].size()); i<limit; ++i)
				{
					depthData<<vector_vs[mpos->second][i]<<delim;
					if (i && i%length == 0) depthData<<'\n';
				}
			}
		}
	}
	else if(depthsingle)
	{
		length = DEPTH_SINGLE_LENGTH_PER_ROW;
		temp = dPos;
		ofstream depthData(temp.c_str());
		for(mpos = chrmaps.begin(); mpos != chrmaps.end(); ++mpos)
		{
			++ptcount;
			depthData<<'>'<<mpos->first<<"\n";
			if(onlycover)
			{
				for(register unsigned long i(1),limit(vector_vs[mpos->second].size()); i<limit; ++i)
				{
					if(vector_vs[mpos->second][i])
						depthData<<1<<delim;
					else
						depthData<<0<<delim;
					if (i && i%length == 0) depthData<<'\n';
				}
			}
			else
			{
				for(register unsigned long i(1),limit(vector_vs[mpos->second].size()); i<limit; ++i)
				{
					depthData<<vector_vs[mpos->second][i]<<delim;
					if (i && i%length == 0) depthData<<'\n';
				}
			}
			depthData<<endl;
		}
	}
}

void OutputCoverageBin()
{
	progress_display ptcount(chrmaps.size(), cout, "Output Coverage to files (Binary)...\n");
	strmaps::iterator mpos;
	ofstream depthData(dPos.c_str(), ios_base::binary | ios_base::out | ios_base::trunc);
	ulong len_tmp(vector_vs.size());
	//Output the quantity of segments
	depthData.write((char*)(&len_tmp), sizeof(ulong));


	for(mpos = chrmaps.begin(); mpos != chrmaps.end(); ++mpos)
	{
		++ptcount;
		//Output the length of the name
		len_tmp = mpos->first.size();
		depthData.write((char*)(&len_tmp), sizeof(ulong));
		//Output the name
		depthData.write(mpos->first.c_str(), mpos->first.size());
		//Output the length of the segment
		len_tmp = (vector_vs[mpos->second].size()-1);
		depthData.write((char*)(&len_tmp), sizeof(ulong));
		if(onlycover)
		{
			ushort s1(1);
			ushort s0(0);
			for(register unsigned long i(1),limit(vector_vs[mpos->second].size()); i<limit; ++i)
				//Output the single point detail
			{
				if(vector_vs[mpos->second][i])
					depthData.write((char*)(&s1), sizeof(ushort));
				else
					depthData.write((char*)(&s0), sizeof(ushort));
			}
		}
		else
		{
			for(register unsigned long i(1),limit(vector_vs[mpos->second].size()); i<limit; ++i)
				//Output the single point detail
				depthData.write((char*)(&vector_vs[mpos->second][i]), sizeof(ushort));
		}
	}
}

void CalcCoverage()
{
	progress_display ptcount(chrmaps.size(), cout, "Calculating Coverage...\n");
	//OutData<<"Coverage of each:\n";
	strmaps::iterator mpos;
	for(mpos = chrmaps.begin(); mpos != chrmaps.end(); ++mpos)
	{
		OutData<<mpos->first<<": ";
		int singleBP(0);
		int effectiveBP(0);
		ulong depth_of_seq(0);
		for(register int i(1),limit(vector_vs[mpos->second].size()); i<limit; ++i)
		{
			if(i<trim && i>=(limit-trim))
				continue;
			if(vector_vs[mpos->second][i]==65535)
				continue;
			++effectiveBP;
			if(vector_vs[mpos->second][i] < minimumDepth)
				continue;
			if(vector_vs[mpos->second][i] > maximumDepth)
				continue;
			if(vector_vs[mpos->second][i])
			{
				depth_of_seq += vector_vs[mpos->second][i];
				++singleBP;
			}
			++distribution[vector_vs[mpos->second][i]];
		}
		countBP += singleBP;
		OutData<<singleBP<<'/'<<effectiveBP<<"\t"<<setprecision(precisionOffset)<<((double)singleBP/effectiveBP*100)<<"\t"<<setprecision(precisionOffset)<<depth_of_seq/(double)singleBP<<endl;
		++ptcount;
	}
	OutData<<endl;
}

void SP()
{
	cout<<"Marking S/P ratio area..."<<endl;
	//cout<<"Int Max: "<<intmax<<endl;
	strmaps::iterator mpos;
	ifstream SPIN(spIn.c_str());
	if(!SPIN)
	{
		cerr<<"Error opening S/P ratio input file: "<<spIn<<endl;
		exit(EXIT_FAILURE);
	}
	string stemp;
	int mark(1);
	if(!pesupport) sp_sui_single.push_back(make_pair<string, uint>("", 0));
	sp_sui_soap.push_back(make_pair<string, uint>("", 0));
	for(;;)
	{
		assert(mark == sp_sui_soap.size());
		string stemp1;
		uint itemp2(0),itemp3(0);
		SPIN>>stemp1>>itemp2>>itemp3;
		if(!SPIN) break;
		mpos = chrmaps.find(stemp1);
		if(mpos == chrmaps.end())
		{
			cerr<<stemp1<<" from S/P ratio input not found in reference!"<<endl;
			exit(EXIT_FAILURE);
		}
		if(!(mark == sp_sui_single.size() && sp_sui_single.size() == sp_sui_soap.size()) && !pesupport)
		{
			cerr<<"!(mark == sp_sui_single.size() && sp_sui_single.size() == sp_sui_soap.size() && !pesupport)"<<endl;
			exit(EXIT_FAILURE);
		}
		if(!pesupport) sp_sui_single.push_back(make_pair<string, uint>(stemp1 + "\t" + lexical_cast<string>(itemp2) + "\t" + lexical_cast<string>(itemp3), 0));
		sp_sui_soap.push_back(make_pair<string, uint>(stemp1 + "\t" + lexical_cast<string>(itemp2) + "\t" + lexical_cast<string>(itemp3), 0));
		if(itemp3 >= sp_vvs[mpos->second].size())
			itemp3 = sp_vvs[mpos->second].size() - 1;
		for(register uint j(itemp2); j<=itemp3; ++j)
			sp_vvs[mpos->second][j] = mark;
		++mark;
	}
	cerr<<"Size of S/P ratio segments: "<<sp_sui_soap.size()<<endl;
	cerr<<"Size of S/P ratio scaffolds: "<<sp_vvs.size()<<endl;
}

void SPOutput()
{
	ofstream O(spOut.c_str());
	if(!O)
	{
		cerr<<"Error opening S/P ratio output file!"<<endl;
		exit(EXIT_FAILURE);
	}
	for(register size_t i(1); i < sp_sui_single.size(); ++i)
		O<<sp_sui_single[i].first<<'\t'<<sp_sui_single[i].second<<'\t'<<sp_sui_soap[i].second<<'\t'<<(sp_sui_soap[i].second?(double)sp_sui_single[i].second/sp_sui_soap[i].second:0)<<endl;
}

void PESupportOutput()
{
	ofstream O(spOut.c_str());
	if(!O)
	{
		cerr<<"Error opening SPE support output file!"<<endl;
		exit(EXIT_FAILURE);
	}
	for(register size_t i(1); i < sp_sui_soap.size(); ++i)
		O<<sp_sui_soap[i].second<<'\t'<<sp_sui_soap[i].first<<endl;
}

void CalcCDS()
{
	cout<<"Summarizing CDS Coverage..."<<endl;
	strmaps::iterator mpos;
	ifstream CDS(cdsPos.c_str());
	if(!CDS)
	{
		cerr<<"Error opening cds input file: "<<cdsPos<<endl;
		exit(EXIT_FAILURE);
	}
	ofstream CDSO;
	if(cdsplot)
	{
		CDSO.open(cdsPos_out.c_str());
		if(!CDSO)
		{
			cerr<<"Error opening cds output file: "<<cdsPos_out<<endl;
			exit(EXIT_FAILURE);
		}
	}
	ofstream CDS_DETAIL;
	if(cdsdetail)
		CDS_DETAIL.open(cdsdetailPos.c_str());
	if(cdsdetail && !CDS_DETAIL)
	{
		cerr<<"Error opening cds details output file: "<<cdsdetailPos<<endl;
		exit(EXIT_FAILURE);
	}
	string stemp;
	vector<int> cds_distribution(65536,0);
	int cds_total(0);
	int cds_covered(0);
	string stemp1,stemp2,stemp3,name;
	for(;;)
	{
		CDS>>stemp1>>stemp2>>stemp3;
		if(!CDS) break;
		name = stemp1 + "\t" + stemp2 + "\t" + stemp3;
		int itemp2(0),itemp3(0);
		itemp2 = atoi(stemp2.c_str());
		itemp3 = atoi(stemp3.c_str());
		mpos = chrmaps.find(stemp1);
		if(mpos == chrmaps.end()) continue;
		stringstream ss;
		//int count(0);
		int cvg_total(0);
		int cvs_covered_sub(0);
		for(register uint j(itemp2); j<=itemp3; ++j)
		{
			if(j < 1 || j>= vector_vs[mpos->second].size())
			{
				cerr<<stemp1<<'\t'<<stemp2<<'\t'<<stemp3<<'\t'<<" Overflow! Omitted."<<endl;
				break;
			}
			if(vector_vs[mpos->second][j] != 65535 && vector_vs[mpos->second][j] > minimumDepth && vector_vs[mpos->second][j] < maximumDepth)
			{
				cvg_total+=vector_vs[mpos->second][j];
				++cds_total;
				if(vector_vs[mpos->second][j])
				{
					++cvs_covered_sub;
					++cds_covered;
				}
				//vector_vs[mpos->second][j] = 65535;
			}
			if(cdsdetail) ss<<vector_vs[mpos->second][j]<<' ';
			++cds_distribution[vector_vs[mpos->second][j]];
		}
		if(cdsdetail)
		{
			CDS_DETAIL<<">"<<name<<'\t'<<itemp3-itemp2+1<<'\t'<<cvs_covered_sub<<'\t'<<setprecision(6)<<((double)cvg_total/cvs_covered_sub)<<endl;
			CDS_DETAIL<<ss.str()<<endl;
		}
	}

	cout<<"Coverage of specific regions: "<<(double)cds_covered/cds_total*100<<"%"<<endl;
	OutData<<"Coverage of specific regions: "<<(double)cds_covered/cds_total*100<<"%"<<endl;
	if(cdsplot)
		for(register int i(cds_plot_lower) ; i<=cds_plot_upper; ++i)
			if(cds_distribution[i])
				CDSO<<setw(5)<<i<<setw(12)<<cds_distribution[i]<<endl;
}

void CalcPoint()
{
	cout<<"Summarizing Point Coverage..."<<endl;
	strmaps::iterator mpos;
	ifstream CDS(pointPos.c_str());
	if(!CDS)
	{
		cerr<<"Error opening point input file: "<<pointPos<<endl;
		exit(EXIT_FAILURE);
	}
	ofstream CDSO(pointOut.c_str());
	if(!CDSO)
	{
		cerr<<"Error opening point output file: "<<pointOut<<endl;
		exit(EXIT_FAILURE);
	}

	string stemp;
	string stemp1,stemp2,name;
	for(;;)
	{
		CDS>>stemp1>>stemp2;
		if(!CDS) break;
		int j(0);
		j = atoi(stemp2.c_str());
		mpos = chrmaps.find(stemp1);
		if(mpos == chrmaps.end() || j < 1 || j>= vector_vs[mpos->second].size())
		{
			CDSO << stemp1 << "\t" << stemp2 << "\t-1\n";
			continue;
		}
		CDSO << stemp1 << "\t" << stemp2 << "\t" << vector_vs[mpos->second][j] << "\n";
	}
}


void AddCvg()
{
	strmaps::iterator mpos;
	ifstream CVG(cPos.c_str());
	if(!CVG)
	{
		cerr<<"Error opening Coverage file: "<<cPos<<endl;
		exit(EXIT_FAILURE);
	}
	//Judge what kind of file (text or binary)
	char firstChar(0);
	CVG.read(&firstChar, 1);
	if(firstChar == '>')
	{
		CVG.seekg(0);
		progress_display ptcount(chrmaps.size(), cout, "Importing Coverage (Text)...\n");
		string cdstemp;
		for(register unsigned long i(1);;++i)
		{
			CVG>>cdstemp;
			if(!CVG) break;
			if(cdstemp[0] == '>')
			{
				++ptcount;
				i = 0;
				cdstemp = cdstemp.substr(1, string::npos);
				mpos = chrmaps.find(cdstemp);
				if(mpos == chrmaps.end())
				{
					cerr<<"Sequence "<<cdstemp<<" missing!"<<endl;
					cerr<<"Wrong input depth file or reference file?"<<endl;
					exit(EXIT_FAILURE);
				}
				continue;
			}
			else
				vector_vs[mpos->second][i] = atoi(cdstemp.c_str());
		}
	}
	else
	{
		CVG.close();
		CVG.clear();
		CVG.open(cPos.c_str(), ios_base::in | ios_base::binary);
		ulong chr_quantity(0);
		ulong len_tmp(0);
		CVG.read((char*)(&chr_quantity), sizeof(ulong));
		if(chr_quantity != chrmaps.size())
		{
			cerr<<"Segment quantity in reference file is different from depthinput file!"<<endl
				<<" Reference: "<<chrmaps.size()<<endl
				<<"Depthinput: "<<chr_quantity<<endl;
			exit(EXIT_FAILURE);
		}
		progress_display ptcount(chrmaps.size(), cout, "Importing Coverage (Binary)...\n");
		string cdstemp;
		for(int x(0); x < chr_quantity; ++x,++ptcount)
		{
			CVG.read((char*)(&len_tmp), sizeof(ulong));
			char* name_tmp = new char[len_tmp+1];
			name_tmp[len_tmp] = '\0';
			CVG.read(name_tmp, len_tmp);
			cdstemp = name_tmp;
			delete[] name_tmp;
			mpos = chrmaps.find(cdstemp);
			if(mpos == chrmaps.end())
			{
				cerr<<"Sequence "<<cdstemp<<" missing!"<<endl;
				cerr<<"Wrong input depth file or reference file?"<<endl;
				exit(EXIT_FAILURE);
			}
			CVG.read((char*)(&len_tmp), sizeof(ulong));
			++len_tmp;
			for(register int i(1); i < len_tmp; ++i)
				CVG.read((char*)(&vector_vs[mpos->second][i]), sizeof(ushort));
		}
	}
}

void AddCvgFromSamtools()
{
	strmaps::iterator mpos;
	igzstream CVG(cPos.c_str());
	if(!CVG)
	{
		cerr<<"Error opening Samtools Coverage file: "<<cPos<<endl;
		exit(EXIT_FAILURE);
	}
	//Judge what kind of file (text or binary)
	progress_display ptcount(chrmaps.size(), cout, "Importing Coverage (Samtools)...\n");
	string prevChr;
	string chr; size_t pos, depth;
	while(true)
	{
		CVG >> chr >> pos >> depth;
		if(!CVG) break;
		if(chr != prevChr)
		{
			prevChr = chr;
			++ptcount;
			mpos = chrmaps.find(chr);
			if(mpos == chrmaps.end())
			{
				cerr<<"Sequence "<<chr<<" missing!"<<endl;
				cerr<<"Wrong input depth file or reference file?"<<endl;
				exit(EXIT_FAILURE);
			}
		}
		if(depth >= 65534)
			depth = 65534;
		vector_vs[mpos->second][pos] = depth;
	}
}

void AddN()
{
	ifstream N(nPos.c_str());
	if(!N)
	{
		cerr<<"Error opening Gap file: "<<nPos<<endl;
		exit(EXIT_FAILURE);
	}
	cerr<<"Masking N gaps..."<<endl;
	string stemp1; int itemp2(0), itemp3(0);
	strmaps::iterator mpos;
	for(;;)
	{
		N>>stemp1>>itemp2>>itemp3;
		if(!N) break;
		mpos = chrmaps.find(stemp1);
		if(mpos == chrmaps.end())
		{
			cerr<<"Error founding "<<stemp1<<endl;
			continue;
		}
		for(register uint j(itemp2); j < itemp3; ++j)
		{
			--totalBP;
			vector_vs[mpos->second][j]=65535;
		}
	}
}

void AddCover()
{
	cout<<"Summarizing Coverage..."<<endl;
	CThreadManager threadMgr(NUM_THREADS);
    threadMgr.SetTimer(0.5);
	vector<param_struct> vparam(iPos.size());
	for(register uint i(0); i<iPos.size(); ++i)
	{
		vparam[i].pos = iPos[i];
		if(plain)
			threadMgr.AddThread( _AddCover_plain, (void*)&vparam[i]);
		else if(general)
			threadMgr.AddThread( _AddCover_general, (void*)&vparam[i]);
		else if(sam)
			threadMgr.AddThread( _AddCover_sam, (void*)&vparam[i]);
		else if(pslsub)
			threadMgr.AddThread( _AddCover_pslsub, (void*)&vparam[i]);
		else if(pslquery)
			threadMgr.AddThread( _AddCover_pslquery, (void*)&vparam[i]);
		else if(maq)
			threadMgr.AddThread( _AddCover_maq, (void*)&vparam[i]);
		else if(m8subject)
			threadMgr.AddThread( _AddCover_m8subject, (void*)&vparam[i]);
		else if(m8query)
			threadMgr.AddThread( _AddCover_m8query, (void*)&vparam[i]);
		else if(mummerquery)
			threadMgr.AddThread( _AddCover_mummerquery, (void*)&vparam[i]);
		else if(axtoitg)
			threadMgr.AddThread( _AddCover_axtoitg, (void*)&vparam[i]);
		else if(axtoiq)
			threadMgr.AddThread( _AddCover_axtoiq, (void*)&vparam[i]);
		else if(qmap)
			threadMgr.AddThread( _AddCover_qmap, (void*)&vparam[i]);
		else
			threadMgr.AddThread( _AddCover_Sub, (void*)&vparam[i]);
	}
	threadMgr.Run();
}

void* _AddCover_plain(void* param)
{
	string ctemp, ctemp2;
	string temp;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());
	while(InData)
	{
		unsigned long end(0),start(0);
		for(;InData;)
		{
			InData>>temp>>ctemp2>>ctemp;
			mpos = chrmaps.find(temp); semapos = chr_sema.find(temp);
			if(mpos == chrmaps.end())
				continue;
			start = atol(ctemp2.c_str());
			end = atol(ctemp.c_str());
			if(end < start)
				swap(start, end);
			if(tag)
				end = start;
			break;
		}
		if(!InData) break;
		pthread_mutex_lock(semapos->second);
		if(end > vector_vs[mpos->second].size())
			end = vector_vs[mpos->second].size();
		for(; start < end ; ++start)
		{
			if(vector_vs[mpos->second][start]<65534)
				++(vector_vs[mpos->second][start]);
		}
		pthread_mutex_unlock(semapos->second);
	}
}

void* _AddCover_general(void* param)
{
	uint rname(idCol - 1), pos(startCol - 1), seq(endCol - 1);
	uint maxCol = max(rname, pos);
	maxCol = max(maxCol, seq);
	string str;
	vector<string> vec_str;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());
	while(InData)
	{
		uint end(0),start(0);
		for(;InData;)
		{
			getline(InData, str);
			if(!InData)
				break;
			if(!str.empty() && str[0]=='#')
				continue;
			vec_str.clear();
			regex_split(back_inserter(vec_str), str);
			if(vec_str.size() <= maxCol)
				continue;
			mpos = chrmaps.find(vec_str[rname]); semapos = chr_sema.find(vec_str[rname]);
			if(mpos == chrmaps.end())
				continue;
			start = atoi(vec_str[pos].c_str());
			(generalLen) ? (end = start + atoi(vec_str[seq].c_str())) : (end = atoi(vec_str[seq].c_str()));
			if(tag)
				end = start;
			break;
		}
		if(!InData) break;
		pthread_mutex_lock(semapos->second);
		if(end > vector_vs[mpos->second].size())
			end = vector_vs[mpos->second].size();
		for(; start < end ; ++start)
		{
			if(vector_vs[mpos->second][start]<65534)
				++(vector_vs[mpos->second][start]);
		}
		pthread_mutex_unlock(semapos->second);
	}
}


#define SINGLE_UNMAPED	0x0004
#define PAIR_UNMAPED	0x0008

void* _AddCover_sam(void* param)
{
	enum str_col{qname=0, flag, rname, pos, mapq, cigar, mrnm, mpos2, isize, seq, qual, tagmark, vtype};
	string str;
	vector<string> vec_str;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());
	while(InData)
	{
		uint end(0),start(0);
		for(;InData;)
		{
			getline(InData, str);
			if(!InData)
				break;
			if(!str.empty() && str[0]=='@')
				continue;
			vec_str.clear();
			regex_split(back_inserter(vec_str), str);
			mpos = chrmaps.find(vec_str[rname]); semapos = chr_sema.find(vec_str[rname]);
			if(atoi(vec_str[flag].c_str()) & (SINGLE_UNMAPED ))//| PAIR_UNMAPED))
				continue;
			if(mpos == chrmaps.end())
				continue;
			start = atoi(vec_str[pos].c_str());
			end = start + vec_str[seq].size() - 1;
			if(tag)
				end = start;
			break;
		}
		if(!InData) break;
		pthread_mutex_lock(semapos->second);
		if(end > vector_vs[mpos->second].size())
			end = vector_vs[mpos->second].size();
		for(; start < end ; ++start)
		{
			if(vector_vs[mpos->second][start]<65534)
				++(vector_vs[mpos->second][start]);
		}
		pthread_mutex_unlock(semapos->second);
	}
}

void* _AddCover_Sub(void* param)
{
	string ctemp, ctemp2, uniqmark, mismatchmark;
	string temp;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;

	imethod InData(pParam->pos.c_str());

	vpsu * sp_sui = &sp_sui_single;
	if(sp)
	{
		for(register size_t i(0); i < isPos.size(); ++i)
			if(isPos[i] == pParam->pos)
				sp_sui = &sp_sui_single;
		for(register size_t i(0); i < ipPos.size(); ++i)
			if(ipPos[i] == pParam->pos)
				sp_sui = &sp_sui_soap;
		if(sp_sui == NULL)
			cerr<<(pParam->pos)<<" not found in both il_single and il_soap!"<<endl;
	}

	vector<uint> mismatch_bypass_vector;
	while(true)
	{
		uint length(0),start(0), mismatch(0);
		for(;;)
		{
			//4 or uniqmark
			InData>>ctemp;InData>>ctemp;InData>>ctemp;InData>>uniqmark;
			if(!InData) break;
			if(onlyuniq && uniqmark!="1")
			{
				InData.ignore(100000,'\n');
				continue;
			}
			//6 for length
			InData>>ctemp;InData>>ctemp2;
			//8 for chr
			InData>>ctemp; InData>>temp;
			mpos = chrmaps.find(temp); semapos = chr_sema.find(temp);
			if(mpos == chrmaps.end())
			{
				InData.ignore(100000,'\n');
				if(!nowarning)
					cerr<<"Omitted one row cuz \""<<temp<<"\" not found!"<<endl;
				continue;
			}
			//9 for start
			InData>>ctemp;
			if(ctemp.size() > 9)
			{
				InData.ignore(100000,'\n');
				if(!nowarning)
					cerr<<"Omitted one row cuz position\""<<ctemp<<"\" make no sense!"<<endl;
				continue;	//Validation
			}
			//10 for mismatchmark
			InData>>mismatchmark;
			start = atoi(ctemp.c_str());
			length = atoi(ctemp2.c_str());
			//bool continue_match(false);
			if(precise)
			{
				mismatch = atoi(mismatchmark.c_str());
				mismatch_bypass_vector.clear();
				if((mismatch == 1) || (mismatch == 2))
				{
					string mismatchpos;
					getline(InData, mismatchpos, '\n');
					vector<string> vec_mismatchpos;
					regex_split(back_inserter(vec_mismatchpos), mismatchpos);
					//mismatchpos = vec_mismatchpos[vec_mismatchpos.size()-1];

					mismatch = SeperateCigar(vec_mismatchpos[vec_mismatchpos.size()-1], mismatch_bypass_vector);
					if(!mismatch)
					{
						if(!nowarning)
							cerr<<"Regex match error on string: "<<endl
								<<temp<<'\t'<<start<<'\t'<<length<<'\t'<<mismatchpos<<endl;
						mismatch_bypass_vector.clear();
						break;
					}
				}
				else
					InData.ignore(100000,'\n');
			}
			else
				InData.ignore(100000,'\n');
			break;
		}
		if(!InData) break;
		if(!length)
		{
			if(!nowarning)
				cerr<<"Omitted one row cuz length equals to 0!"<<endl;
			continue;				//Validation
		}
		if(tag)
			length = 1;
		pthread_mutex_lock(semapos->second);

		uint j(start+length);
		if(j > vector_vs[mpos->second].size())
			j = vector_vs[mpos->second].size();
		for(uint loopcount(0); start<j ; ++start,++loopcount)
		{
			if(precise)
				for(uint checkmismatch_loopcount(0); checkmismatch_loopcount < mismatch_bypass_vector.size(); ++checkmismatch_loopcount)
					if(mismatch_bypass_vector[checkmismatch_loopcount] == loopcount)
						continue;
			if(vector_vs[mpos->second][start]<65534)
				++(vector_vs[mpos->second][start]);
		}

		if(sp)
		{
			uint end(start+length-1);
			if(end >= sp_vvs[mpos->second].size())
				end = sp_vvs[mpos->second].size() - 1;
			//if(sp_vvs[mpos->second][start] >= (*sp_sui).size() || sp_vvs[mpos->second][end] >= (*sp_sui).size())
			//{
				//cerr<<sp_vvs[mpos->second][start]<<'\t'<<sp_vvs[mpos->second][end]<<'\t'<<((*sp_sui).size())<<endl;
			//	continue;
			//}
			if(sp_vvs[mpos->second][start] == 0 && sp_vvs[mpos->second][end] == 0)
				continue;
			else if(sp_vvs[mpos->second][start] != 0 && sp_vvs[mpos->second][end] == 0)
				++((*sp_sui)[sp_vvs[mpos->second][start]]).second;
			else if(sp_vvs[mpos->second][start] == 0 && sp_vvs[mpos->second][end] != 0)
				++((*sp_sui)[sp_vvs[mpos->second][end]]).second;
			else if(sp_vvs[mpos->second][start] == sp_vvs[mpos->second][end] && sp_vvs[mpos->second][start] != 0)
				++((*sp_sui)[sp_vvs[mpos->second][start]]).second;
			else if(sp_vvs[mpos->second][start] != sp_vvs[mpos->second][end] && sp_vvs[mpos->second][start] != 0 && sp_vvs[mpos->second][end] != 0)
			{
				++((*sp_sui)[sp_vvs[mpos->second][start]]).second;
				++((*sp_sui)[sp_vvs[mpos->second][end]]).second;
			}
			else
			{
				cerr<<sp_vvs[mpos->second][start]<<'\t'<<sp_vvs[mpos->second][end]<<" Exceptional S/P ratio condition!"<<endl;
			}
		}
		pthread_mutex_unlock(semapos->second);
	}
}

void* _AddCover_pslsub(void* param)
{
	string ctemp, ctemp2;
	string length, substart;
	vector<string> length_vec, substart_vec;
	regex e(",");
	string temp;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());
	while(InData)
	{
		uint end(0),start(0);
		for(;InData;)
		{
			//Omitting data
			InData>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp;
			//14 for chr
			InData>>temp;
			mpos = chrmaps.find(temp); semapos = chr_sema.find(temp);
			if(mpos == chrmaps.end())
			{
				InData.ignore(100000,'\n');
				continue;
			}
			InData>>ctemp>>ctemp>>ctemp>>ctemp>>length;
			InData>>ctemp>>substart;
			length_vec.clear();
			substart_vec.clear();
			regex_split(back_inserter(length_vec), length, e);
			regex_split(back_inserter(substart_vec), substart, e);
			break;
		}
		if(!InData) break;
		pthread_mutex_lock(semapos->second);
		for(uint vec_it(0); vec_it < length_vec.size(); ++vec_it)
		{
			start = atoi(substart_vec[vec_it].c_str());
			end = start + atoi(length_vec[vec_it].c_str());
			if(end >= vector_vs[mpos->second].size())
				end = vector_vs[mpos->second].size();
			for(; start < end ; ++start)
			{
				if(vector_vs[mpos->second][start]<65534)
					++(vector_vs[mpos->second][start]);
			}
		}
		pthread_mutex_unlock(semapos->second);
	}
}

void* _AddCover_pslquery(void* param)
{
	string ctemp, ctemp2;
	string length, substart;
	vector<string> length_vec, substart_vec;
	regex e(",");
	string temp;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());
	while(InData)
	{
		uint end(0),start(0);
		for(;InData;)
		{
			//Omitting data
			InData>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp;
			//14 for chr
			InData>>temp;
			mpos = chrmaps.find(temp); semapos = chr_sema.find(temp);
			if(mpos == chrmaps.end())
			{
				InData.ignore(100000,'\n');
				continue;
			}
			InData>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>length;
			InData>>substart>>ctemp;
			length_vec.clear();
			substart_vec.clear();
			regex_split(back_inserter(length_vec), length, e);
			regex_split(back_inserter(substart_vec), substart, e);
			break;
		}
		if(!InData) break;
		pthread_mutex_lock(semapos->second);
		for(uint vec_it(0); vec_it < length_vec.size(); ++vec_it)
		{
			start = atoi(substart_vec[vec_it].c_str());
			end = start + atoi(length_vec[vec_it].c_str());
			if(end >= vector_vs[mpos->second].size())
				end = vector_vs[mpos->second].size();
			for(; start < end ; ++start)
			{
				if(vector_vs[mpos->second][start]<65534)
					++(vector_vs[mpos->second][start]);
			}
		}
		pthread_mutex_unlock(semapos->second);
	}
}

void* _AddCover_maq(void* param)
{
	string ctemp, ctemp2;
	string temp;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());
	while(InData)
	{
		uint end(0),start(0);
		for(;InData;)
		{
			//Omitting data
			InData>>ctemp;
			//2 for chr
			InData>>temp;
			mpos = chrmaps.find(temp); semapos = chr_sema.find(temp);
			if(mpos == chrmaps.end())
			{
				InData.ignore(100000,'\n');
				continue;
			}
			//3 for start
			InData>>ctemp2;
			//14 for end
			InData>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp;
			start = atoi(ctemp2.c_str());
			end = atoi(ctemp.c_str()) + start - 1;
			InData.ignore(100000,'\n');
			if(tag)
				end = start;
			break;
		}
		if(!InData) break;
		pthread_mutex_lock(semapos->second);
		for(register uint limit(vector_vs[mpos->second].size()); start <= end ; ++start)
		{
			if(start >= limit) break;
			if(vector_vs[mpos->second][start]<65534)
				++(vector_vs[mpos->second][start]);
		}
		pthread_mutex_unlock(semapos->second);
	}
}

void* _AddCover_m8subject(void* param)
{
	string ctemp, ctemp2;
	string temp;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());
	while(InData)
	{
		uint end(0),start(0);
		for(;InData;)
		{
			//Omitting data
			InData>>ctemp;
			//2 for chr
			InData>>temp;
			mpos = chrmaps.find(temp); semapos = chr_sema.find(temp);
			if(mpos == chrmaps.end())
			{
				InData.ignore(100000,'\n');
				continue;
			}
			//9 for start
			InData>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp2;
			//10 for end
			InData>>ctemp;
			start = atoi(ctemp2.c_str());
			end = atoi(ctemp.c_str());
			if(end < start)
				swap(start, end);
			InData.ignore(100000,'\n');
			break;
		}
		if(!InData) break;
		pthread_mutex_lock(semapos->second);
		for(register uint limit(vector_vs[mpos->second].size()); start <= end ; ++start)
		{
			if(start >= limit) break;
			if(vector_vs[mpos->second][start]<65534)
				++(vector_vs[mpos->second][start]);
		}
		pthread_mutex_unlock(semapos->second);
	}
}

void* _AddCover_m8query(void* param)
{
	string ctemp, ctemp2;
	string temp;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());
	while(InData)
	{
		uint end(0),start(0);
		for(;InData;)
		{
			//1 for chr
			InData>>temp;
			mpos = chrmaps.find(temp); semapos = chr_sema.find(temp);
			if(mpos == chrmaps.end())
			{
				InData.ignore(100000,'\n');
				continue;
			}
			//7 for start
			InData>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp>>ctemp2;
			//8 for end
			InData>>ctemp;
			start = atoi(ctemp2.c_str());
			end = atoi(ctemp.c_str());
			if(start > end)
				swap(start, end);
			InData.ignore(100000,'\n');
			break;
		}
		if(!InData) break;
		pthread_mutex_lock(semapos->second);
		for(register uint limit(vector_vs[mpos->second].size()); start <= end ; ++start)
		{
			if(start >= limit) break;
			if(vector_vs[mpos->second][start]<65534)
				++(vector_vs[mpos->second][start]);
		}
		pthread_mutex_unlock(semapos->second);
	}
}

void* _AddCover_mummerquery(void* param)
{
	string ctemp, ctemp2;
	string temp(">");
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());

	bool four_column(false);
	//mummer format judgement
	for(;;)
	{
		while(temp.find(">") != string::npos)
			getline(InData, temp);
		stringstream ss(temp);
		int i(0);
		for(;;++i)
		{
			ss>>temp;
			if(!ss)
				break;
		}
		if(i == 4)
			four_column = true;
		InData.seekg(0);
		break;
	}

	//Data processing
	while(InData)
	{
		uint end(0),start(0);
		bool reverse_mark(false);
		bool read_signal(true);
		string name;
		for(;;)
		{
			//For special mark or ref or ref position
			if(read_signal)
				InData>>temp;
			else
				read_signal = true;

			if(!InData) break;
			if(!temp.empty() && temp[0] == '>')
			{
				reverse_mark = false;
				InData>>temp;
				name = temp;
				mpos = chrmaps.find(name); semapos = chr_sema.find(name);
				if(mpos == chrmaps.end())
				{
					cerr<<"Fatal error! Name "<<name<<" in mummer result file missing in reference, input file bypassed!"<<endl;
					return (void*)(NULL);
				}
				InData>>temp;
				if(temp.find("Reverse") != string::npos)
				{
					reverse_mark = true;
					InData>>temp;
				}
				if(temp[0] == '>')
				{
					read_signal = false;
					continue;
				}
			}

			if(four_column)
			{
				//For query position
				InData>>ctemp;
				InData>>ctemp;
			}
			else
			//For query position
				InData>>ctemp;

			//For length
			InData>>ctemp2;

			if(atoi(ctemp2.c_str()) < mummer_limit)
				continue;
			/*if(atoi(ctemp2.c_str()) > 100000)
			{
				cerr<<"Record "<<ctemp<<'\t'<<ctemp2<<" > 100000"<<"! Omited."<<endl;
				continue;
			}*/

			if(reverse_mark)
			{
				end = vector_vs[mpos->second].size() - atoi(ctemp.c_str());
				start = end - atoi(ctemp2.c_str());
			}
			else
			{
				start = atoi(ctemp.c_str());
				end = start + atoi(ctemp2.c_str());
			}
			break;
		}

		if(tag)
			end = start + 1;

		if(!InData) break;
		pthread_mutex_lock(semapos->second);
		for(register uint limit(vector_vs[mpos->second].size()); start < end ; ++start)
		{
			if(start >= limit) break;
			if(vector_vs[mpos->second][start]<65534)
				++(vector_vs[mpos->second][start]);
		}
		pthread_mutex_unlock(semapos->second);
	}
}

void* _AddCover_qmap(void* param)
{
	string temp(">"), strTemp;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());
	ulong position(0);
	uint met(0), uMet(0);
	bool locked(false);

	//Data processing
	while(InData)
	{
		getline(InData, temp, '\n');
		if(!InData)
			break;
		if(!temp.empty() && temp[0] == '>')
		{
			if(locked)
				pthread_mutex_unlock(semapos->second);
			locked = false;
			position = 0;
			string name(temp.substr(1, temp.find_first_of(" \t") - 1));
			mpos = chrmaps.find(name); semapos = chr_sema.find(name);
			if(mpos == chrmaps.end())
			{
				cerr<<"Fatal error! Name "<<name<<" in qmap missing in reference, input file bypassed!"<<endl;
				return (void*)(NULL);
			}
			pthread_mutex_lock(semapos->second);
			locked = true;
		}
		else
		{
			stringstream ss(temp);
			ss>>strTemp>>met>>uMet;
			vector_vs[mpos->second][position] += (met + uMet);
		}
	}
}

void* _AddCover_axtoitg(void* param)
{
	string ctemp, ctemp2;
	string temp;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());
	while(InData)
	{
		uint end(0),start(0);
		for(;InData;)
		{
			//2 for chr
			InData>>temp;
			if(InData && (temp.empty() || isalpha(temp[0])))
				continue;
			InData>>temp;
			mpos = chrmaps.find(temp); semapos = chr_sema.find(temp);
			if(mpos == chrmaps.end())
			{
				InData.ignore(100000,'\n');
				continue;
			}
			//3 for start
			InData>>ctemp2;
			//4 for end
			InData>>ctemp;
			start = atoi(ctemp2.c_str());
			end = atoi(ctemp.c_str());
			InData.ignore(100000,'\n');
			break;
		}
		if(!InData) break;
		pthread_mutex_lock(semapos->second);
		for(register uint limit(vector_vs[mpos->second].size()); start <= end ; ++start)
		{
			if(start >= limit) break;
			if(vector_vs[mpos->second][start]<65534)
				++(vector_vs[mpos->second][start]);
		}
		pthread_mutex_unlock(semapos->second);
	}
}

void* _AddCover_axtoiq(void* param)
{
	string ctemp, ctemp2;
	string temp;
	string strand;
	strmaps::iterator mpos;
	semaphore_maps::iterator semapos;
	param_struct* pParam = (param_struct*) param;
	static ushort progress(1);
	cout<<progress++<<'/'<<iPos.size()<<": "<<pParam->pos<<endl;
	imethod InData(pParam->pos.c_str());
	while(InData)
	{
		uint end(0),start(0);
		for(;InData;)
		{
			//5 for chr
			InData>>temp;
			if(InData && (temp.empty() || isalpha(temp[0])))
				continue;
			InData>>temp>>temp>>temp>>temp;
			mpos = chrmaps.find(temp); semapos = chr_sema.find(temp);
			if(mpos == chrmaps.end())
			{
				InData.ignore(100000,'\n');
				continue;
			}
			//6 for start
			InData>>ctemp2;
			//7 for end
			InData>>ctemp;
			//8 for strand
			InData>>strand;
			if(!strand.empty() && strand[0]=='-')
			{
				start = vector_vs[mpos->second].size() - atoi(ctemp.c_str());
				end = vector_vs[mpos->second].size() - atoi(ctemp2.c_str());
			}
			else
			{
				start = atoi(ctemp2.c_str());
				end = atoi(ctemp.c_str());
			}
			InData.ignore(100000,'\n');
			break;
		}
		if(!InData) break;
		pthread_mutex_lock(semapos->second);
		for(register uint limit(vector_vs[mpos->second].size()); start <= end ; ++start)
		{
			if(start >= limit) break;
			if(vector_vs[mpos->second][start]<65534)
				++(vector_vs[mpos->second][start]);
		}
		pthread_mutex_unlock(semapos->second);
	}
}

ulong AddTotalBP()
{
	ulong singleBP(0);
	string ctemp;
	ifstream refData(temp.c_str());
	if(!refData)
	{
		cerr<<"Error while opening: "<<temp<<endl;
		exit(EXIT_FAILURE);
	}
	getline(refData,ctemp);		//Abandon the first information line
	while(refData>>ctemp) singleBP+=ctemp.size();
	totalBP += (singleBP - trim * 2);
	return (singleBP+1);
}

ulong AddTotalBP_refsingle()
{
	return
		refsingle_marker[temp];
}

ulong AddTotalBP_fastat()
{
	return
		refsingle_marker[temp];
}

void BuildMemory()
{
	progress_display ptcount(chr.size(), cout, "Building Memory Blocks...\n");
	strsets::iterator pos;
	for(pos=chr.begin(); pos!=chr.end(); ++pos)
	{
		++ptcount;
		ulong singleBP(0);
		if(refsingle)
		{
			temp = *pos;
			singleBP=AddTotalBP_refsingle();
		}
		else if(reffastat)
		{
			temp = *pos;
			singleBP=AddTotalBP_fastat();
		}
		else if(::ref)
		{
			temp = rPos + '/' + *pos + ".fa";
			singleBP=AddTotalBP();
		}

		vector_vs.push_back(vs(singleBP,0));
		if(sp)
			sp_vvs.push_back(vi(singleBP, 0));
	}
}

void MakeMaps()
{
	cout<<"Shadowing Map..."<<endl;
	strsets::iterator pos;
	for(pos=chr.begin(); pos!=chr.end(); ++pos)
	{
		static int i(0);
		chrmaps.insert(make_pair(*pos,i));
		chr_sema.insert(make_pair(*pos, new pthread_mutex_t));
		if(pthread_mutex_init(chr_sema[*pos], NULL))
		{
			cerr<<"Mutex on "<<*pos<<" initialization error!"<<endl;
			exit(EXIT_FAILURE);
		}
		i++;
	}
	cerr<<"Mutex Lock created: "<<chr_sema.size()<<endl;
}

void ParaDeal(int & argc,char** argv)
{
	string temp;
	char ctemp[10000];
	if(argc == 1)
	{
		cerr<<USAGE;
		exit(EXIT_SUCCESS);
	}
	cout<<endl<<"Parameters List:";
	//Parameters Deals
	for(register int i(1),next(0); i<argc; i++)
	{
		int rightPara(0);
		if (next)
		{
			cout<<argv[i]<<' ';
			--next;
			continue;
		}
		if (strcmp(argv[i],"-h") == 0 || strcmp(argv[i],"-help") == 0)
		{
			cerr<<USAGE;
			exit(EXIT_SUCCESS);
		}
		if (strcmp(argv[i], "-precisionOffset") == 0 || strcmp(argv[i], "-po") == 0)
		{
			ComfirmAllDigit(argv[i+1]);
			precisionOffset = atoi(argv[i+1]);
			rightPara++;
			next++;
		}
		if (strcmp(argv[i], "-minimumDepth") == 0 || strcmp(argv[i], "-mind") == 0)
		{
			ComfirmAllDigit(argv[i+1]);
			minimumDepth = atoi(argv[i+1]);
			rightPara++;
			next++;
		}
		if (strcmp(argv[i], "-maximumDepth") == 0 || strcmp(argv[i], "-maxd") == 0)
		{
			ComfirmAllDigit(argv[i+1]);
			maximumDepth = atoi(argv[i+1]);
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-nowarning") == 0 || strcmp(argv[i],"-nw") == 0)
		{
			nowarning = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-onlyuniq") == 0 || strcmp(argv[i],"-ou") == 0)
		{
			onlyuniq = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-precise") == 0)
		{
			precise = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-plain") == 0)
		{
			plain = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-sam") == 0)
		{
			sam = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-qmap") == 0)
		{
			qmap = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-pslquery") == 0)
		{
			pslquery = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-pslsub") == 0)
		{
			pslsub = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-maq") == 0)
		{
			maq = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-m8subject") == 0)
		{
			m8subject = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-m8query") == 0)
		{
			m8query = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-mummerquery") == 0)
		{
			mummerquery = true;
			ComfirmAllDigit(argv[i+1]);
			mummer_limit = atoi(argv[i+1]);
			if(mummer_limit < 0) mummer_limit = 20;
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-axtoitg") == 0)
		{
			axtoitg = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-axtoiq") == 0)
		{
			axtoiq = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-cvg") == 0)
		{
			if(phy)
			{
				cerr<<"Parameter \"-cvg\" can't be used with \"-phy\""<<endl;
				exit(EXIT_FAILURE);
			}
			cvg = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-nospace") == 0)
		{
			nospace = true;
		}
		if (strcmp(argv[i],"-phy") == 0)
		{
			if(cvg)
			{
				cerr<<"Parameter \"-phy\" can't be used with \"-cvg\""<<endl;
				exit(EXIT_FAILURE);
			}
			phy = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-tag") == 0)
		{
			if(phy)
			{
				cerr<<"Parameter \"-tag\" can't be used with \"-phy\""<<endl;
				exit(EXIT_FAILURE);
			}
			cvg = true;
			tag = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-nocalc") == 0 || strcmp(argv[i],"-nc") == 0)
		{
			nocalc = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-onlycover") == 0 || strcmp(argv[i],"-oc") == 0)
		{
			onlycover = true;
			rightPara++;
		}
		if (strcmp(argv[i],"-p") == 0)
		{
			ComfirmAllDigit(argv[i+1]);
			NUM_THREADS = atoi(argv[i+1]);
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-plot") == 0)
		{
			plot = true;
			if(i + 3 >= argc)
			{
				cerr<<"\n-plot parameter format error!\n";
				exit(1);
			}
			plotPos = argv[i+1];
			ComfirmAllDigit(argv[i+2]);
			ComfirmAllDigit(argv[i+3]);
			plot_lower = atoi(argv[i+2]);
			plot_upper = atoi(argv[i+3]);
			if(plot_lower < 0) plot_lower = 0;
			if(plot_upper > 65534) plot_upper = 65534;
			rightPara++;
			next+=3;
		}
		if (strcmp(argv[i],"-general") == 0 || strcmp(argv[i],"-generallen") == 0)
		{
			if(strcmp(argv[i],"-general") == 0)
				general = true;
			else
				general = generalLen = true;
			if(i + 3 >= argc)
			{
				cerr<<"\n-general or -generalLen parameter format error!\n";
				exit(1);
			}
			ComfirmAllDigit(argv[i+1]);
			ComfirmAllDigit(argv[i+2]);
			ComfirmAllDigit(argv[i+3]);
			idCol = atoi(argv[i+1]);
			startCol = atoi(argv[i+2]);
			endCol = atoi(argv[i+3]);
			rightPara++;
			next+=3;
		}
		if (strcmp(argv[i],"-duplicate") == 0)
		{
			ComfirmAllDigit(argv[i+1]);
			duplicate = atoi(argv[i+1]);
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-insertupper") == 0 || strcmp(argv[i],"-iu") == 0)
		{
			ComfirmAllDigit(argv[i+1]);
			insert_Limit = atoi(argv[i+1]);
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-insertlower") == 0 || strcmp(argv[i],"-ilo") == 0)
		{
			ComfirmAllDigit(argv[i+1]);
			insert_LowerLimit = atoi(argv[i+1]);
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-addn") == 0)
		{
			nPos = argv[i+1];
			addN=true;
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-depthinput") == 0 || strcmp(argv[i],"-di") == 0)
		{
			cPos = argv[i+1];
			depthinput=true;
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-depthinputsam") == 0 || strcmp(argv[i],"-dis") == 0)
		{
			cPos = argv[i+1];
			depthinputsam=true;
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-cdsinput") == 0 || strcmp(argv[i],"-ci") == 0)
		{
			cdsPos = argv[i+1];
			cdsinput=true;
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-point") == 0)
		{
			pointPos = argv[i+1];
			pointOut = argv[i+2];
			rightPara++;
			next+=2;
		}
		if (strcmp(argv[i],"-sp") == 0)
		{
			sp = true;
			spIn = argv[i+1];
			spOut = argv[i+2];
			rightPara++;
			next+=2;
		}
		if (strcmp(argv[i],"-pesupport") == 0 || strcmp(argv[i],"-ps") == 0)
		{
			pesupport = true;
			sp = true;
			spIn = argv[i+1];
			spOut = argv[i+2];
			rightPara++;
			next+=2;
		}
		if (strcmp(argv[i],"-cdsplot") == 0 || strcmp(argv[i],"-cp") == 0)
		{
			cdsplot = true;
			cdsPos_out = argv[i+1];
			ComfirmAllDigit(argv[i+2]);
			ComfirmAllDigit(argv[i+3]);
			cds_plot_lower = atoi(argv[i+2]);
			cds_plot_upper = atoi(argv[i+3]);
			if(cds_plot_lower < 0) cds_plot_lower = 0;
			if(cds_plot_upper > 65534) cds_plot_upper = 65534;
			rightPara++;
			next+=3;
		}
		if (strcmp(argv[i],"-cdsdetail") == 0 || strcmp(argv[i],"-cd") == 0)
		{
			cdsdetail =  true;
			cdsdetailPos = argv[i+1];
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-window") == 0)
		{
			window = true;
			windowPos_out = argv[i+1];
			ComfirmAllDigit(argv[i+2]);
			window_size = atoi(argv[i+2]);
			if(window_size < 0) window_size = 0;
			rightPara++;
			next+=2;
		}
		if (strcmp(argv[i],"-trim") == 0)
		{
			ComfirmAllDigit(argv[i+1]);
			trim = atoi(argv[i+1]);
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-depth") == 0)
		{
			dPos = argv[i+1];
			depth=true;
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-depthsingle") == 0 || strcmp(argv[i],"-ds") == 0)
		{
			dPos = argv[i+1];
			depthsingle=true;
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-depthsinglebin") == 0 || strcmp(argv[i],"-dsb") == 0)
		{
			dPos = argv[i+1];
			depthsinglebin=true;
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-o") == 0)
		{
			oPos = argv[i+1];
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-i") == 0)
		{
			addon = true;
			for(register int j(1); (i+j)<argc; ++j)
			{
				if (argv[i+j][0] == '-') break;
				igzstream InData_sub(argv[i+j]);
				if(!InData_sub)
				{
					cerr<<"Error opening Soap file:"<<argv[i+j]<<endl;
					exit(EXIT_FAILURE);
				}
				iPos.push_back(argv[i+j]);
				next++;
			}
			rightPara++;
		}
		if (strcmp(argv[i],"-il") == 0)
		{
			addon = true;
			igzstream InData(argv[i+1]);
			if(!InData)
			{
				cerr<<"Error opening Soap list-file:"<<argv[i+1]<<endl;
				exit(EXIT_FAILURE);
			}
			for(InData>>ctemp; InData; InData>>ctemp)
			{
				iPos.push_back(ctemp);
			}
			for(register uint i(0); i<iPos.size(); ++i)
			{
				igzstream InData_sub(iPos[i].c_str());
				if(!InData_sub)
				{
					cerr<<"Error opening Soap file:"<<iPos[i]<<endl;
					exit(EXIT_FAILURE);
				}
			}
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-il_single") == 0)
		{
			addon = true;
			il_single = true;
			igzstream InData(argv[i+1]);
			if(!InData)
			{
				cerr<<"Error opening Soap single list-file:"<<argv[i+1]<<endl;
				exit(EXIT_FAILURE);
			}
			for(InData>>ctemp; InData; InData>>ctemp)
			{
				iPos.push_back(ctemp);
				isPos.push_back(ctemp);
			}
			for(register uint i(0); i<isPos.size(); ++i)
			{
				igzstream InData_sub(isPos[i].c_str());
				if(!InData_sub)
				{
					cerr<<"Error opening Soap single file:"<<isPos[i]<<endl;
					exit(EXIT_FAILURE);
				}
			}
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-il_soap") == 0)
		{
			addon = true;
			il_soap = true;
			igzstream InData(argv[i+1]);
			if(!InData)
			{
				cerr<<"Error opening Soap PE list-file:"<<argv[i+1]<<endl;
				exit(EXIT_FAILURE);
			}
			for(InData>>ctemp; InData; InData>>ctemp)
			{
				iPos.push_back(ctemp);
				ipPos.push_back(ctemp);
			}
			for(register uint i(0); i<ipPos.size(); ++i)
			{
				igzstream InData_sub(ipPos[i].c_str());
				if(!InData_sub)
				{
					cerr<<"Error opening Soap PE file:"<<ipPos[i]<<endl;
					exit(EXIT_FAILURE);
				}
			}
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-ref") == 0)
		{
			::ref = true;
			igzstream InData(argv[i+1]);
			if(!InData)
			{
				cerr<<"Error opening Reference list file: "<<argv[i+1]<<endl;
				exit(EXIT_FAILURE);
			}
			for(InData>>temp; InData; InData>>temp)
			{
				rPos = temp.substr(0,temp.rfind('/'));
				temp=temp.substr((temp.rfind('/')+1), string::npos);
				temp=temp.substr(0, temp.find(".fa"));
				chr.insert(temp);
			}
			rightPara++;
			next++;
		}
		if (strcmp(argv[i],"-refsingle") == 0 || strcmp(argv[i],"-rs") == 0)
		{
			rPos = argv[i+1];
			rightPara++;
			next++;
			refsingle=true;
		}
		if (strcmp(argv[i],"-reffastat") == 0 || strcmp(argv[i],"-rfs") == 0)
		{
			rPos = argv[i+1];
			rightPara++;
			next++;
			reffastat=true;
		}
		if (rightPara == 0)
		{
			cout<<"Parameter error on \""<<argv[i]<<"\""<<endl;
			exit(EXIT_FAILURE);
		}
		else
			cout<<endl<<argv[i]<<' ';
	}
	cout<<endl<<"# End of parameters list"<<endl<<endl;

	//Parameter Logicals
	if((cvg && phy) || (!cvg && !phy))
	{
		cerr<<"At least and only one of the \"-cvg\", \"-phy\" or \"-tag\" should be selected!"<<endl;
		exit(EXIT_FAILURE);
	}
	if(pslsub && phy)
	{
		cerr<<"-pslsub should only be used with -cvg other than -phy!"<<endl;
		exit(EXIT_FAILURE);
	}
	if(pslquery && phy)
	{
		cerr<<"-pslquery should only be used with -cvg other than -phy!"<<endl;
		exit(EXIT_FAILURE);
	}
	if(maq && phy)
	{
		cerr<<"-maq should only be used with -cvg other than -phy!"<<endl;
		exit(EXIT_FAILURE);
	}
	if(maq && phy)
	{
		cerr<<"-m8 should only be used with -cvg other than -phy!"<<endl;
		exit(EXIT_FAILURE);
	}
	if(depth && depthsingle)
	{
		cerr<<"Only one of the \"-depth\" or \"-depthsingle\" could be selected!"<<endl;
		exit(EXIT_FAILURE);
	}
	if(depthsinglebin && depthsingle)
	{
		cerr<<"Only one of the \"-depthsinglebin\" or \"-depthsingle\" could be selected!"<<endl;
		exit(EXIT_FAILURE);
	}
	/*if(oPos.empty())
	{
		cerr<<"\"-o\" should be defined!"<<endl;
		exit(EXIT_FAILURE);
	}*/
	if(duplicate && depthinput)
	{
		cerr<<"\"-duplication\" can't be used with \"-depthinput\"!"<<endl;
		exit(EXIT_FAILURE);
	}
	if(sp && !pesupport && (!il_single || !il_soap))
	{
		cerr<<"\"-sp\" should be used with \"-il_single\" & \"il_soap\"!"<<endl;
		exit(EXIT_FAILURE);
	}
	if(sp && phy && !pesupport)
	{
		cerr<<"\"-sp\" should not be used with \"-phy\""<<endl;
		exit(EXIT_FAILURE);
	}
	if(!phy && pesupport)
	{
		cerr<<"\"-pesupport\" should be used with \"-phy\""<<endl;
		exit(EXIT_FAILURE);
	}
	if(refsingle)
		AddRefSingle();
	if(reffastat)
		AddFastat();

	//Opening General Output file
	if(oPos.empty())
		cerr<<"General output (-o) undefined, disabled!"<<endl;
	else
	{
		OutData.open(oPos.c_str());
		if(!OutData)
		{
			cerr<<"Error opening output file!"<<endl;
			exit(EXIT_FAILURE);
		}
	}
}

void AddRefSingle()
{
	string temp;
	string ctemp;
	igzstream InData(rPos.c_str());
	if(!InData)
	{
		cerr<<"Error opening Reference file!";
		exit(1);
	}
	cout<<"Picking out segments from reference file..."<<endl
		<<"Number of segments:"<<setw(10)<<'0'<<flush;
	for(register uint i(0); InData ;++i)
	{
		ulong seg_length(0);
		if(!i)
			getline(InData, ctemp, '\n');
		if(ctemp[0] == '>')
		{
			temp=ctemp;
			temp=temp.substr(1,(temp.find_first_of(" \t")-1));
			while(true)
			{
				getline(InData, ctemp, '\n');
				if(!InData) break;
				if(ctemp[0] == '>') break;
				seg_length += ctemp.size();
			}
			totalBP += (seg_length - trim * 2);
			if(!nowarning && refsingle_marker.find(temp) != refsingle_marker.end())
			{
				cerr<<"Error! There are repeating segment name in reference!"<<endl;
				exit(EXIT_FAILURE);
			}
			refsingle_marker.insert(make_pair<string, ulong>(temp, seg_length + 1));
			cout<<"\b\b\b\b\b\b\b\b\b\b"
				<<setw(10)<<refsingle_marker.size()<<flush;
		}
		else
		{
			cerr<<"Marker \">\" not found in line "<<ctemp<<"in refsingle!"<<endl;
			exit(EXIT_FAILURE);
		}
		chr.insert(temp);
	}
	InData.clear();
	InData.close();
	cout<<endl;
}

void AddFastat()
{
	string temp;
	string ctemp;
	igzstream InData(rPos.c_str());
	if(!InData)
	{
		cerr<<"Error opening Reference file!";
		exit(1);
	}
	cout<<"Picking out segments from FaStat file..."<<endl
		<<"Number of segments:"<<setw(10)<<'0'<<flush;
	for(register uint i(0); InData ;++i)
	{
		string name, size, size_withoutn;
		for(;;)
		{
			InData>>name>>size>>size_withoutn;
			if(!InData)
				break;
			//name=name.substr(1,(temp.find_first_of(" \t")-1));
			totalBP += (atol(size.c_str()) - trim * 2);
			if(!nowarning && refsingle_marker.find(temp) != refsingle_marker.end())
			{
				cerr<<"Error! There are repeating segment name in FaStat!"<<endl;
				exit(EXIT_FAILURE);
			}
			refsingle_marker.insert(make_pair<string, ulong>(name, atol(size.c_str()) + 1));
			cout<<"\b\b\b\b\b\b\b\b\b\b"
				<<setw(10)<<refsingle_marker.size()<<flush;
			chr.insert(name);
		}
	}
	InData.clear();
	InData.close();
	cout<<endl;
}

void Output_list()
{
	int i(1);
	cout<<"Output List:"<<endl;
	cout<<setw(2)<<i++<<": "<<oPos<<endl;
	if(depth)
		cout<<setw(2)<<i++<<": "<<dPos<<"/*"<<endl;
	if(depthsingle)
		cout<<setw(2)<<i++<<": "<<dPos<<endl;
	if(cdsinput)
		cout<<setw(2)<<i++<<": "<<cdsPos_out<<endl;
}

void ComfirmAllDigit(const char* rhs)
{
	for(register uint i(0); i<strlen(rhs); ++i)
	{
		if(!isdigit(rhs[i]))
		{
			cerr<<rhs<<" is not a number!"<<endl;
			exit(EXIT_FAILURE);
		}
	}
}

inline uint SeperateCigar(string & str, vector<uint> & vec)
{
	string pos;
	uint character(0);
	uint accumulate(0);
	for(register int i(0); i < str.size(); ++i)
	{
		if(isalpha(str[i]) || str[i] == '-')
		{
			if(!pos.empty())
				accumulate += atoi(pos.c_str());
			vec.push_back(accumulate + character);
			pos.clear();
			++character;
		}
		else
			pos += str[i];
	}

	return character;
}

