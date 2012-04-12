/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-8-9
 *CLASS NAME: FileListManager
 *FUNCTION : a class that can control file list. support open, readline, close etc
 *			file functions.
 *FILE NAME : FileListManager.h
 *UPDATE DATE : 2010-8-19
 *UPDATE BY : Bill Tang
 *UPDATA DATA : 2010-8.30
 *UPDATA BY :ZHOU GUYUE
 *UPDATE: 2010-11-2 change the m_soap_file_vec's type, to vector<*igzstream>.
					change the interface which is related to soap file.
 *UPDATA BY :BIll Tang
 *******************************************************************************
 */

#pragma once

/*define some error number*/
#define CNSFILE_ERROR		-2		//can not open cnsFile
#define BASEFILE_ERROR		-1		//can not open listFile
#define SOAPFILE_ERROR		-2		//can not open soapFile
#define SAMFILE_ERROR		-2		//can not open sam/bamFile 
#define SETLIMIT_ERROR		-3		//can not set limit

#define SETLIMIT_OK          3
#define INDEX_OVER_FLOW		NULL	// index is over flow
#define COUNT_ERROR			-3		// the readwin_vec size not equal to the out files number.
#define INPUT_ERROR			-4		// the input files has something wrong.
#define FILE_END			2		// file is end.

#include "SamCtrl.h"
#include "Readwin.h"
#include "soap_snp.h"
#include "CThreadPool.h"
#include "Call_winManager.h"

using namespace std;

class FileListManager
{
	int m_file_number;     // file number in file list.
	my_ifstream m_base_file;	// the handle of the file that contain file list.
	vector<igzstream*>  m_soap_file_vec;	// the handle vector of the files which are soap format.
	vector<SamCtrl*>  m_sam_file_vec;	// the handle vector of the files which are sam/bam format.
	vector<gzoutstream*>  m_cns_file_vec;	// the handle vector of the output files.
	vector<fstream*> m_matrix_file_vec;	// a vector stored the file handle which are matrix files
	int m_stop_pos;		// record the position once one of the samples reach a break point that the chromosome changed.
	int m_endfile_count;				// record the file's count that reach the end.
public:	
	FileListManager();
	virtual ~FileListManager();
	// open output file function
	virtual int openCnsFile(const string& path, const string& outdir, const char* mode, Parameter * para);
	// open all the file in the file list
	virtual int openAli(const string& listfilename, const char* mode);
	// open sam/bam file
	virtual int openSamFile(const string& listfilename, const char* mode);
	// open soap format file
	virtual int openSoapFile(const string& listfilename);
	// read a line from the *vec[index]
	virtual int readLine(string& line, const int index);
	// close alignment files
	virtual void closeAliFiles(void);
	// close consensus files
	virtual void closeCnsFiles(void);
	// read lines to the readwin_vec. The lines' start position must inside a win size. When the lines' end position exceed the win,put it to the next win.
	virtual int readWin(vector<Readwin>& readwin_vec, Parameter * para);
	// get the soap file handle pointer m_soap_file_vec[index]
	virtual igzstream * getSoapFile(const int index);
	// get the soap file handle pointer m_soap_file_vec[index]
	virtual SamCtrl* getSamFile(const int index);
	// read record from soapfile and send to the readwin.
	virtual int readWinFromSoap(Readwin& readwin, igzstream& soapfile);
	// read record from samCtrl and send to the readwin.
	virtual int readWinFromSam(Readwin& readwin, SamCtrl& samCtrl);
	// get the file number
	virtual int getFileNum(void);
	// get handle of consensus file.
	virtual gzoutstream* getCnsFile(const int index);
	// open matrix files with open_mode mode
	virtual int openMatrixFile(const std::string matrix_list, std::ios_base::open_mode mode, const std::string alignment_list = "");
	// get the matrix file handle
	virtual fstream* getMatrixFile(const int index);
	//set limit of file number;
	//virtual int setlimit(int filenumber); //update by zhukai on 2010-12-28
	
};

//2010-11-08

typedef struct _read_win_args
{
	//FileListManager * fileListManager_p;
	Readwin * readwin_p;
	igzstream * m_soap_file_p;
	SamCtrl * m_sam_file_p;
	int * ret_p;
	int index;

	inline _read_win_args(Readwin * b, igzstream *c, SamCtrl * f,int * d, int e)
	{
		//fileListManager_p = a;
		readwin_p = b;
		m_soap_file_p = c;
		m_sam_file_p = f;
		//cerr <<"111 ret_p = d "<<endl;
		ret_p = d;
		index = e; 
	};
	/*inline _read_win_args(Readwin * b,SamCtrl *c, int * d, int e)
	{
		//fileListManager_p = a;
		readwin_p = b;
		m_sam_file_p = c;
		m_soap_file_p = NULL;
		ret_p = d;
		index = e;
	};*/
	inline _read_win_args()
	{
		//fileListManager_p = NULL;
		readwin_p = NULL;
		m_soap_file_p = NULL;
		m_sam_file_p = NULL;
		ret_p = NULL;
		index = 0;
	};
}READ_WIN_ARGS;

class Read_win_Task : public CTask 
{
public:
	FileListManager fileListManager;
	// overload the function Run().
	int Run();
	/*{
		_call_winManager_soap2cns_norm(this->m_ptrData);
		return 1;
	}*/
};

// 2010-11-12  the function readWin's parameter structure.
typedef struct _big_read_win_args
{
	FileListManager * fileListManager;
	vector<Readwin> * readwin_vec;
	Parameter * para;

	inline _big_read_win_args(FileListManager * a, vector<Readwin> * b, Parameter * c)
	{
		fileListManager = a;
		readwin_vec = b;
		para = c;
	};

	inline _big_read_win_args()
	{
		fileListManager = NULL;
		readwin_vec = NULL;
		para = NULL;
	};
} BIG_READ_WIN_ARGS;

void *_flieListManager_readWin(void * args);
