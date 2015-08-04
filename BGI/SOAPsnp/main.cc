
#include "soap_snp.h"
#include <getopt.h>
#include <ctime>
#include "FileListManager.h"
#include "MatrixManager.h"
#include "Readwin.h"
#include "Call_winManager.h"
#include "accessControl.h"
#include "SfsMethod.h"
using namespace std;
#include <istream>
#include <sys/mman.h>

#define VERSIONSTR "v2.0"
 
/**
 * DATE:  
 * FUNCTION:  return the help information
 * PARAMETER:  void
 * RETURN:  help information
 */

int usage() {

	cerr<<"		****************************************************************			"<<endl;
	cerr<<"		*		  Software name :SOAPsnp                       *			"<<endl;
	cerr<<"		*	     	        Version : v 2.0                        *			"<<endl;
	cerr<<"		*		    Last Update : 2010.10.15		       *			"<<endl;
	cerr<<"		*		         Author : Bill Tang                    *			"<<endl;
	cerr<<"		*		         E_mail : tangzhoubiao@genomics.org.cn *			"<<endl;
	cerr<<"		*				  zhouguyue@genomics.org.cn    *			"<<endl;
	cerr<<"		*				  zhukai@genomics.org.cn       *			"<<endl;
	cerr<<"		*		      Copyright : BGI. All Rights Reserved.    *			"<<endl;
	cerr<<"		****************************************************************			"<<endl;
	cerr<<"`Compulsory Parameters:"<<endl;
	cerr<<"		-i <FILE> Input SORTED Soap Result(Can be .gz file)"<<endl;
	cerr<<"		-S <FILE> Input SORTED SAM Result"<<endl;
	cerr<<"		-B <FILE> Input SORTED BAM Result"<<endl;
	cerr<<"		-l if the Soap is filelist input"<<endl;
	cerr<<"		-d <FILE> Reference Sequence in fasta format"<<endl;
	cerr<<"		-o <FILE> (<DIR> if input is file list, if want to call sfs, you don't need to set this option, but you want call both, you can set this option) Output consensus file Optional Parameters:(Default in [])"<<endl;
	cerr<<"		-z <Char> ASCII chracter standing for quality==0 [@]"<<endl;
	cerr<<"		-g <Double> Global Error Dependency Coefficient, 0.0(complete dependent)~1.0(complete independent)[0.9]"<<endl;
	cerr<<"		-p <Double> PCR Error Dependency Coefficient, 0.0(complete  dependent)~1.0(complete independent)[0.5]"<<endl;
	cerr<<"		-r <Double> novel altHOM prior probability [0.0005]"<<endl;
	cerr<<"		-e <Double> novel HET prior probability [0.0010]"<<endl;
	cerr<<"		-t set transition/transversion ratio to 2:1 in prior probability"<<endl;
	cerr<<"		-s <FILE> Pre-formated dbSNP information"<<endl;
	cerr<<"		-2 specify this option will REFINE SNPs using dbSNPs information [Off]"<<endl;
	cerr<<"		-a <Double> Validated HET prior, if no allele frequency known [0.1]"<<endl;
	cerr<<"		-b <Double> Validated altHOM prior, if no allele frequency known[0.05]"<<endl;
	cerr<<"		-j <Double> Unvalidated HET prior, if no allele frequency known [0.02]"<<endl;
	cerr<<"		-k <Double> Unvalidated altHOM rate, if no allele frequency known[0.01]"<<endl;
	cerr<<"		-u Enable rank sum test to give HET further penalty for better accuracy.[Off]"<<endl;
	//cerr<<"-n Enable binomial probability calculation to give HET for better accuracy. [Off]"<<endl;
	cerr<<"		-m Enable monoploid calling mode, this will ensure all consensus as HOM  and you probably should SPECIFY"<<endl<<
	      "		   higher altHOM rate. [Off]"<<endl;	
	cerr<<"		-q Only output potential SNPs. Useful in Text output mode. [Off]"<<endl;
	cerr<<"		-M <FILE>(<DIR> if set -l parameter) Output the quality calibration matrix;the matrix can be reused with -I "<<endl;
	cerr<<"	           if you rerun the program"<<endl;
	cerr<<"		-I <FILE>(<FILELIST> if set -l parameter) Input previous quality calibration matrix. It cannot be used "<<endl<<
	      "		   simutaneously with -M ."<<endl;
	cerr<<"		-L <short> maximum length of read [45] ."<<endl;
	cerr<<"		-Q <short> maximum FASTQ quality score [40] ."<<endl;
//	cerr<<"		-F <int> Output format. 0: Text; 1: GLFv2; 2: GPFv2.[0]"<<endl;
	cerr<<"		-E <String> Extra headers EXCEPT CHROMOSOME FIELD specified in GLFv2 output. "<<endl<<
	      "	           Format is \"TypeName1:DataName1:TypeName2:DataName2\"[""] ."<<endl;
	cerr<<"		-T <FILE> Only call consensus on regions specified in FILE. Format: ChrName\\tStart\\tEnd."<<endl;
	cerr <<"		-C <int> The CPU number which you want to set(It only work with -l).[4]" << endl;
	cerr <<"		-Z <int> Define if the output consensus is .gz format. 1 for .gz format, 0 for non .gz format.[0]" << endl;
	cerr << "		-f <FILE> the sfs parameter file. The file format: (must be used with -l parameter)" << endl << endl;
	cerr << "			doBay	0/1	" << endl;
	cerr << "			doJoint	0/1		" << endl;
	cerr << "			writeFr 0/1			" << endl;
	cerr << "			outfiles	<FILE PATH>	" << endl;
	cerr << "			underFlowProtect	0/1	" << endl;
	cerr << "			allowPhatZero	0/1	" << endl;
	cerr << "			bias	[0,1]	" << endl;
	cerr << "			pvar	[0,1]	" << endl;
	cerr << "			eps	[0,1]	" << endl;
	cerr << "			pThres 	[0,1]	" << endl;
	cerr << "			sigleMinorJoint	0/1" << endl;
	cerr << "			sigleMinorBay	0/1" << endl ;
	cerr <<"			qs		<int>quality score threshold" << endl ; 
	cerr <<"			lightoutput    0/1" << endl <<endl; 

	//cerr<<"-S <FILE> Output summary of consensus"<<endl;
	cerr << "		-h Display this help" << endl;
	exit(1);
	return 0;
}

/**
 * DATE:  
 * FUNCTION:  get help information function
 * PARAMETER:  void
 * RETURN:  help information
 */
int readme() {
	return usage();
}

/**
 * DATE:  
 * FUNCTION:  get the parament from script
 * PARAMETER:  argc :the number of parmament
			   argv[] : the parament 
			   files : point to the files
 * RETURN:   
 */
void get_c ( int argc, char * argv[], Parameter *para, Files *files)
{
	while(( para->c = getopt(argc,argv,"i:d:o:z:g:p:r:e:ts:2a:b:j:k:unmqM:I:L:Q:S:B:F:E:T:h:f:C:Z:l")) != -1) {
		switch(para->c) {
			case 'l':
			{
				/*if file_list is true means alignment is a file list*/
				para ->file_list = true;
				break;
			}	
			case 'i':
			{
				// Soap Alignment Result
				/*files->soap_result.clear();
				//files->soap_result.open(optarg);
				//if( ! files->soap_result) {
				cerr<<"No such file or directory:"<<optarg<<endl;
				exit(1);
				}*/
				//soap alignment path
				para->alignment_name = optarg;
				break;
			}
			case 'd':
			{
				// The reference genome in fasta format
				files->ref_seq.clear();
				files->ref_seq.open(optarg);
				if( ! files->ref_seq) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					exit(1);
				}
				files->ref_seq.clear();
				break;
			}
			case 'o':
			{
				//output file path
				para->consensus_name = optarg;
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
				files->dbsnp.clear();
				files->dbsnp.open(optarg);
				if( ! files->ref_seq) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					exit(1);
				}
				files->dbsnp.clear();
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
				para->matrix_file = optarg;
				
				break;
			}
			case 'I':
			{
				para->matrix_file = optarg;
				
				para->is_matrix_in = true;
				break;
			}
			case 'S':
			{
				//files->summary.open(optarg);
				//// Output the summary of consensus
				//if( ! files->summary ) {
				//	cerr<<"No such file or directory: "<<optarg<<endl;
				//	exit(1);
				//}
				para->alignment_name = optarg;
				para->in_mode[0] = 'r';  // the input file is sam format
				break;
			}
			case 'B':
			{
				para->alignment_name = optarg;
				para->in_mode[0] = 'r';
				para->in_mode[1] = 'b';             //the input file is bam format 
				break;
			}
			case 'L':
			{
				para->read_length = atoi(optarg);     //the read length
				break;
			}
			case 'Q':
			{
				para->q_max = optarg[0];            //set q_max
				if(para->q_max < para->q_min) {
					cerr<< "FASTQ quality character error: Q_MAX > Q_MIN" <<endl;
				}
				break;
			}
			case 'F': {
				para->glf_format = atoi(optarg);        //output glf format 
				break;
			}
			case 'E': {
				para->glf_header = optarg;
				break;
			}
			case 'T': {
				files->region.clear();
				files->region.open(optarg);
				files->region.clear();
				para->region_only = true;          //region_only means have target region file
				break;
			}
			case 'f': {
				para->sfs_path = optarg;            //the path of sfs file   
				break;
			}
			case 'C': 
			{
				para->CPU = atoi(optarg);            //the number of the CPU for thread
				break;
			}
			case 'Z':
			{
				if (atoi(optarg))
				{// set the output cns format
					para->is_cns_gz = true;
				}
				break;
			}
			case 'h':readme();break;
			case '?':usage();break;
			default: cerr<<"Unknown error in command line parameters"<<endl;
		}
	}
	 
 
	if (para->alignment_name == "")
	{
		usage();
		exit(1);
	}
	// change the order between file_list and sfs_path. 2010-12-14
	if (para->file_list == true)
	{
		if (para->sfs_path == "" )
		{
			cerr << "ERROR: -l must be used with -f" << endl;
			exit(0);
		}
	}
	
	if (para->CPU < 1)
	{
		cerr << "\tERROR: You have set the wrong CPU number " << para->CPU << endl
			 << "\tUse -h for more information." << endl;
		exit(0);
	}
}

/**
 * DATE:  
 * FUNCTION:  open soap/bam/sam file and if it's single open matrix file for read or write
 * PARAMETER:  para , files, filelistmanager for file list , sfsmethod
 * RETURN:   
 */
void open_file ( Parameter *para,Files *files, FileListManager *fileListManager, SfsMethod *sfsMethod)
{
	if (para->file_list == true)         //open all the file in the list                    
	{
		para->ret = fileListManager->openAli(para->alignment_name,para->in_mode);   //openAli function return file number or error information
		if (BASEFILE_ERROR == para->ret)
		{
			cerr << "ERROR: File list " << para->alignment_name << "cannot open" << endl;
			exit(1);
		}
		else if (SOAPFILE_ERROR == para->ret)
		{
			cerr << "ERROR: FILE OPEN" << endl;
			exit(1);
		}

		// if need to call sfs, then init the para->sfs
		if (para->sfs_path != "")
		{
			para->sfs = sfs_init();           //new one sfs object
			int ret = init_sfs_para(para->sfs_path, para->sfs);    //init_sfs_para function get the sfs information an return success or error information
			if (ret == SFS_PARA_ERROR)
			{
				exit(1);
			}
			// initialization the SfsMethod :likes , bays, binomLookup matrix, return 0
			sfsMethod->init(fileListManager->getFileNum());

			//open sfs file include frq, bay, sfs . return the open file success or error information
			ret = files->OpenSfsfile(para->sfs->outfiles, para->sfs->writeFr, para->sfs->doBay, para->sfs->doJoint);
			if (ret == -1)
			{
				cerr << "ERROR: FILE OPEN" << para->sfs->outfiles << endl;
				exit(0);
			}
		}

	}//end of file list open 
	else       
	{	
		/*if it is sam\bam open sam_result*/
		if(para->in_mode[0] == 'r')
		{
			files->sam_result.open(para->alignment_name.c_str(), para->in_mode);
		}
		/*open a soap file */
		else 
		{
			files->soap_result.clear();
			files->soap_result.open(para->alignment_name.c_str());
			if( !files->soap_result) 
			{
				cerr<<"No such file or directory:"<<para->alignment_name.c_str()<<endl;
				exit(1);
			}
		}

		if (para->matrix_file != "")
		{
			if (para->is_matrix_in)  
			{
				// Input the calibration matrix , open matrix_file for read 
				files->matrix_file.open(para->matrix_file.c_str(), fstream::in);
				if( ! files->matrix_file) 
				{
					cerr<<"No such file or directory:"<<para->matrix_file<<endl;
					exit(1);
				}
				files->matrix_file.clear();
			}
			else 
			{
				// Output the calibration matrix , open matrix_file for write 
				files->matrix_file.open(para->matrix_file.c_str(), fstream::out);
				if( ! files->matrix_file) 
				{
					cerr<<"Cannot creat file :"<<para->matrix_file<<endl;
					exit(1);
				}
				files->matrix_file.clear();
			}
		}
	}
}

/**
 * DATE:  
 * FUNCTION:  open output format
 * PARAMETER:  para , files, filelistmanager for file list , sfsmethod
 * RETURN:   
 */
void output_format (Parameter *para, FileListManager *fileListManager, Files *files)
{
	if ( ! para->glf_format ) 
	{   //not out put glf format
		if(para->file_list == true)
		{
			para->ret = fileListManager->openCnsFile(para->alignment_name, para->consensus_name, "", para);
			if (CNSFILE_ERROR == para->ret)
			{
				cerr << "ERROR: cannot create consensus file" << endl;
				exit(1);
			}
		}
		else
		{
			// Normal SOAPsnp tab-delimited text format
			//update by zhukai on 2010-12-09
			if(para->is_cns_gz)
			{
				files->consensus = new myogzstream();
				para->consensus_name += ".gz";
			}
			else
			{
				files->consensus = new myofstream();
			}

			files->consensus->clear();			
			files->consensus->open(para->consensus_name.c_str(),std::ios::out);
			if( ! files->consensus->is_open() ) {
				cerr<<"Cannot creat file:" <<para->consensus_name<<endl;
				exit(1);
			}
			files->consensus->clear();
			// set the output format. 2010-12-15 Bill
			files->consensus->set_f(std::ios::showpoint);
		}
	}
	else 
	{
		// SOAPsnp-defined GLF and baseinfo format
		//update by zhukai on 2010-12-09
		if(para->is_cns_gz)
		{
			files->consensus = new myogzstream();
			para->consensus_name += ".gz";
		}
		else
		{
			files->consensus = new myofstream();
		}
		files->consensus->clear();
		files->consensus->open(para->consensus_name.c_str(), ios::binary);
		if(!files->consensus->is_open()) 
		{
			cerr<<"Cannot creat file:" <<para->consensus_name<<endl;
			exit(1);
		}
		files->consensus->clear();

		/*ouput consensus_name.baseinfo*/
		files->baseinfo.clear();
		string baseinfo_name = para->consensus_name + ".baseinfo";
		files->baseinfo.open(baseinfo_name.c_str());
		if(!files->baseinfo) 
		{
			cerr<<"Cannot creat file:" <<baseinfo_name<<endl;
			exit(1);
		}
		files->baseinfo.clear();

		/*ouput consensus_name.index*/
		files->o_region.clear();
		string o_region_name = para->consensus_name + ".index";
		files->o_region.open(o_region_name.c_str());
		if(!files->o_region) 
		{
			cerr<<"Cannot creat file:" <<o_region_name<<endl;
			exit(1);
		}
		files->o_region.clear();
	}
}

/**
 * DATE:  
 * FUNCTION:  set region information , close referece and dbsnp files  
 * PARAMETER:  para , files, filelistmanager for file list , sfsmethod
 * RETURN:   
 */
void read_chr (Genome * genome, Files *files, Parameter *para)
{
	files -> ref_seq.close();
	files -> dbsnp.close();
	clog << "Reading Chromosome and dbSNP information Done."<<endl;
	if (para->region_only && files->region) 
	{
		genome->read_region(files->region, para);
		clog<<"Read target region done."<<endl;
	}
}

/**
 * DATE:  
 * FUNCTION:  deal file list 
 * PARAMETER:  para , files, filelistmanager for file list , sfsmethod
 * RETURN:   
 * UPDATE: 2010-12-14 reduce the judge code if (para->sfs != NULL).
 */
void deal_list(Genome * genome, Parameter *para, FileListManager *fileListManager, SfsMethod *sfsMethod, Files *files)
{
	time_t timep;
	int file_num = fileListManager->getFileNum();

	MatrixManager matrixManager;
	Call_winManager call_winManager(para->read_length, file_num, genome); //initilize call_winManager object
	vector<Readwin> readwin_vec;
	ThreadManager threadManager(para->CPU);         //initilize thread

	//begin to get matrix
	if (para->is_matrix_in)              //matrix is existed
	{
		// open matrix files to input.
		fileListManager->openMatrixFile(para->matrix_file, ios::in);
		//cerr <<"490" <<endl;
		for ( int i=0;i<file_num;++i)
		//for (int i = 0; i < fileListManager->getFileNum(); ++i)   //update by zhukai on 2010-12-09
		{
			if (CREATE_MAT_FAILED == matrixManager.addMatrix(*(fileListManager->getMatrixFile(i)), para))
			{
				cerr << "ERROR: add matrix failed!" << endl;
				exit(0);
			}
			Readwin new_readwin;
			readwin_vec.push_back(new_readwin);
		}
	}
	else  
	{
		matrixManager.setMatrixNum(file_num);
		if (para->matrix_file != "")
		{
			// open matrix file to output.
			fileListManager->openMatrixFile(para->matrix_file, ios::out, para->alignment_name);
		}

		if (para->in_mode[0] == 'r')      //genorate sam/bam file matrix
		{
			for ( int i=0;i<file_num;++i)
		    //for (int i = 0; i < fileListManager->getFileNum(); ++i)   //update by zhukai on 2010-12-09
			{
				//cerr <<"begin sam matrix "<<endl;
				time(&timep);
				clog << ctime(&timep);

				MATRIX_ARGS *tmp_args = new MATRIX_ARGS(&matrixManager, para, genome, 0,
					fileListManager->getMatrixFile(i),
					fileListManager->getSamFile(i),i);
				// add the data to ThreadManager.
				threadManager.AddThread(_matrixManager_addMatrix, tmp_args);

				/*fstream *mat_out_p = fileListManager.getMatrixFile(i);

				if (CREATE_MAT_FAILED == matrixManager.addMatrix(*(fileListManager.getSamFile(i)), para, genome, mat_out_p))
				{
				cerr << "ERROR: add matrix failed!" << endl;
				exit(0);
				}*/

				Readwin new_readwin;
				readwin_vec.push_back(new_readwin);
			}
		}
		else              //genorate soap file matrixs
		{
			for ( int i=0;i<file_num;++i)
			//for (int i = 0; i < fileListManager->getFileNum(); ++i)   //update by zhukai on 2010-12-09
			{
				MATRIX_ARGS *tmp_args = new MATRIX_ARGS(&matrixManager, para, genome,
					fileListManager->getSoapFile(i),
					fileListManager->getMatrixFile(i),
					0,i);
				// add the data to ThreadManager.
				threadManager.AddThread(_matrixManager_addMatrix, tmp_args);

				/*fstream *mat_out_p = fileListManager.getMatrixFile(i);
				if (CREATE_MAT_FAILED == matrixManager.addMatrix(*(fileListManager.getSoapFile(i)), para, genome, mat_out_p))
				{
				cerr << "ERROR: add matrix failed!" << endl;
				exit(0);
				}*/

				Readwin new_readwin;
				readwin_vec.push_back(new_readwin);
			}
		}

		threadManager.Run();
		threadManager.Reset();
	}
	clog << "Correction Matrix Done!" << endl;

	fileListManager->closeAliFiles();
	fileListManager->openAli(para->alignment_name, para->in_mode);

	 
	/***************************************************************************/
	map<Chr_name, Chr_info*>::iterator itr = genome->chromosomes.begin();
	int rw_start = itr->second->getStartPos() - global_win_size; 
	call_winManager.setLast_start(itr->second->getStartPos() );
	for (int i = 0; i < readwin_vec.size(); ++i)
	{
		readwin_vec[i].setWin_pos(rw_start);
	}
	/***************************************************************************/

	int count = 0;
	int sub_count = 0;
	time(&timep);
	clog << "begin to call !!!!!  time : " <<  ctime(&timep);
	CThreadPool threadpool(para->CPU);
	CALL_WIN_ARGS * call_win_args = NULL;
	Call_win_Task * cw_Task;
	// two vector used to release resource.
	vector<CALL_WIN_ARGS*> call_win_args_vec;
	vector<Call_win_Task*> cw_Task_vec;
	pthread_t readwin_pid;
	BIG_READ_WIN_ARGS _read_win_args(fileListManager, &readwin_vec, para);
	pthread_create(&readwin_pid, NULL, _flieListManager_readWin, (void*)&_read_win_args);
	sem_t * sem_read_p = &(para->sem_read);
	sem_t * sem_call_cns_p = &(para->sem_call_cns);
	sem_t * sem_readwin_return_p = &(para->sem_readwin_return);
	int ret;

	//2010-11-16 
	pthread_t call_sfs_pid;
	BIG_CALL_SFS_ARGS _call_sfs_args(sfsMethod, files, para);
	pthread_create(&call_sfs_pid, NULL, _sfsMethod_callsfs, (void*)&_call_sfs_args);
	sem_t * sem_call_sfs_p = &(sfsMethod->sem_call_sfs);
	sem_t * sem_call_sfs_return_p = &(sfsMethod->sem_call_sfs_return);

	do
	{
		//ret = fileListManager->readWin(readwin_vec, para);
		sem_wait(sem_call_cns_p);
		ret = para->ret;		
		if (COUNT_ERROR == ret)
		{
			cerr << "ERROR: something wrong with the file number!" << endl;
			exit(0);
		}
		else if (INPUT_ERROR == ret)
		{
			cerr << "ERROR: something wrong with input files!" << endl;
			exit(0);
		}
		else if (COME_NEW_CHR == ret)
		{
			for ( int i = 0; i < file_num; ++i)
			{
				call_winManager.soap2cns(readwin_vec[i].getReadwin(), fileListManager->getCnsFile(i), files->baseinfo, genome, matrixManager.getMatrix(i), para, i, readwin_vec[i].getSoap2cnsIdx(), *sfsMethod);
				readwin_vec[i].winChange();
				call_winManager.soap2cns(readwin_vec[i].getReadwin(), fileListManager->getCnsFile(i), files->baseinfo, genome, matrixManager.getMatrix(i), para, i, readwin_vec[i].getSoap2cnsIdx(), *sfsMethod);
				readwin_vec[i].reset();
			}
		}

		sem_wait(sem_call_sfs_return_p);

		for ( int i=0;i<file_num;++i)
		{		
			call_win_args = new CALL_WIN_ARGS(&call_winManager, &(readwin_vec[i].getReadwin()), fileListManager->getCnsFile(i), &(files->baseinfo), genome, matrixManager.getMatrix(i), para, i, readwin_vec[i].getSoap2cnsIdx(), sfsMethod);
			call_win_args_vec.push_back(call_win_args);
			cw_Task = new Call_win_Task();
			cw_Task->SetData(call_win_args);
			cw_Task_vec.push_back(cw_Task);
			threadpool.AddTask(cw_Task);
			//call_winManager.soap2cns(readwin_vec[i].getReadwin(), *(fileListManager.getCnsFile(i)), files->baseinfo, genome, matrixManager.getMatrix(i), para, i, readwin_vec[i].getSoap2cnsIdx(), sfsMethod);
		}
		//12-10
		sem_post(sem_readwin_return_p);

		while (threadpool.getTaskSize() != 0)
		{
			usleep(10);
		}

		for (int j = 0; j < cw_Task_vec.size(); ++j)
		{
			delete call_win_args_vec[j];
			delete cw_Task_vec[j];
		}
		call_win_args_vec.clear();
		cw_Task_vec.clear();

		sfsMethod->setidxProcess();
		sem_post(sem_call_sfs_p);
		sem_post(sem_read_p);

	} while(ret != FILE_END);// end of do

	// wait for the last win that call_sfs is processing. 2010-12-14 Bill
	sem_wait(sem_call_sfs_return_p);
	
	for (int i = 0; i < file_num; ++i)
	{
		call_winManager.dealTail(fileListManager->getCnsFile(i), files->baseinfo, genome, matrixManager.getMatrix(i), para, i, *sfsMethod);
	}

	// send last win to call_sfs. and wait for is's end. 2010-12-14 Bill
	sfsMethod->setidxProcess();
	sem_post(sem_call_sfs_p);
	sem_wait(sem_call_sfs_return_p);
	sfsMethod->file_end_flag = 1;
	// make sure that the call_sfs can judge the sfsMethod->file_end_flag and exit. 2010-12-14 Bill
	sem_post(sem_call_sfs_p);
	pthread_join(call_sfs_pid, NULL);

	threadpool.StopAll();
	fileListManager->closeAliFiles();
	fileListManager->closeCnsFiles(); 
	
	time(&timep);
	clog << "Main:" << "end of call!!!" << ctime(&timep);
}

/**
 * DATE:  
 * FUNCTION:  deal single file  
 * PARAMETER:  para , files, filelistmanager for file list , sfsmethod
 * RETURN:   void
 */
void deal_single(Genome * genome,  Files *files, SfsMethod *sfsMethod, Parameter *para)
{
	
	Prob_matrix * mat = new Prob_matrix((para->rank_sum_mode));

	//generate correction matrix 
	if( ! para->is_matrix_in) 
	{     //correction matrix is not exist
		//Read the soap result and give the calibration matrix
		
		if(para->in_mode[0] == 'r')
		{   // generate sam/bam correction matrix
			mat->matrix_gen(files->sam_result, para, genome);
		} else 
		{                        // generate soap correction matrix
			mat->matrix_gen(files->soap_result, para, genome);
		}
		if (files->matrix_file) 
		{
			mat->matrix_write(files->matrix_file, para);
		}
	}
	else 
	{                //correction matrix is exist, just read it 
		mat->matrix_read(files->matrix_file, para);
	}
	files->matrix_file.close();
	clog<<"Correction Matrix Done!"<<endl;
	
	//genorate the prior matrix
	mat->prior_gen(para);
	//do rank sum 
	if (para->rank_sum_mode)
	{
		mat->rank_table_gen();
	}
	
	//intilize call win object with the read_length (-L argument)
	Call_win info(para->read_length);
	info.initialize(0);
	//Call the consensus
	if (para->in_mode[0] != 'r')        //deal with the soap file
	{
		files->soap_result.close();
		files->soap_result.clear();
		files->soap_result.open(para->alignment_name.c_str()); //open the soap output file for output
		files->soap_result.clear();
		info.soap2cns(files->soap_result, files->consensus, files->baseinfo, genome, mat, para, *sfsMethod, 0);
		files->soap_result.close();
		files->consensus->close();
	} 
	else                             //deal with the bam/sam file
	{
		files->sam_result.close();
		files->sam_result.open(para->alignment_name.c_str(), para->in_mode);
		info.soap2cns(files->sam_result, files->consensus, files->baseinfo, genome, mat, para, *sfsMethod, 0);
		files->sam_result.close();
		files->consensus->close();
	}
}

int main ( int argc, char * argv[])
{
	AccessControl("/share/backup/luoruibang/accessControl/SOAPsnp-v2", VERSIONSTR);
	
	Parameter *para = new Parameter();
	FileListManager fileListManager;
	SfsMethod sfsMethod; 
	Files *files = new Files();

	get_c(argc,argv,para,files);
	open_file(para, files, &fileListManager, &sfsMethod);
	
	/*if it is noraml soapsnp format open,else change to binary and open */
	output_format(para, &fileListManager, files);
	time_t timep;
    time(&timep);

    clog << "genome start " << ctime(&timep);
	//Read the chromosomes into memory, get the reference and dbsnp information 
	Genome * genome = new Genome(files->ref_seq, files->dbsnp);
	//read the targ region file information 
	read_chr(genome, files, para);
	time(&timep);
    clog << "genome end " << ctime(&timep);
	
	//begin to deal list or single 
	// change the order between list and single. 2010-12-14
	if (para->file_list != true)
	{
		deal_single(genome, files, &sfsMethod, para);
	}
	else if (para->sfs != NULL)
	{
		time(&timep);
		clog << ctime(&timep);
		deal_list(genome, para, &fileListManager, &sfsMethod, files);
	}
	else
	{
		cerr << "Please check your parameters." << endl;
		exit(0);
	}

	delete files;
	delete para;
	cerr<<"Consensus Done!"<<endl;
	return 0;
}

