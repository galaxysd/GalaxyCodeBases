/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-8-10
 *CLASS NAME: Call_winManager, Call_win_Task
 *FUNCTION : a class that can control Call_win vector.
 *FILE NAME : Call_winManager.h
 *UPDATE DATE : 2010-9-15
 *UPDATE BY : Bill Tang
 *UPDATE: 2010-10-13 add struct CALL_WIN_ARGS and Class Call_win_Task.
 *******************************************************************************
 */

#pragma once
#include "soap_snp.h"
#include "CThreadPool.h"

#define CRE_WIN_FAILED -1	// create Call_win object failed.
#define INDEX_EXCEED	-2	// m_cw_vec index exceed
#define SOAP2CNS_SUCCESS 1	// soap2cns successful.

class Call_winManager
{
	// Call_win objects vector
	vector<Call_win*> m_cw_vec;
	Call_winManager(void);
	
public:
	Call_winManager(ubit64_t read_len, int file_num, Genome* genome);
	virtual ~Call_winManager(void);
	// add a Call_win object to the Call_win vector
	virtual int addCallwin(ubit64_t read_len);
	// soap to consensus
	virtual int soap2cns(vector<Soap_format>& alignment_vec, gzoutstream * consensus, my_ofstream& baseinfo, Genome* genome, Prob_matrix* mat, Parameter* para, int index, int ali_index, SfsMethod &sfsMethod);
	// deal the tail of sample
	virtual int dealTail(gzoutstream * consensus, my_ofstream& baseinfo, Genome* genome, Prob_matrix* mat, Parameter* para, int index, SfsMethod &sfsMethod);
	// set call win's last_start
	void setLast_start(int start_pos);
};

/**
 * DATE: 2010-10-13
 * FUNCTION: a parameter structure that about Call_win.
 * PARAMETER: 
 * RETURN: 
 */
typedef struct _call_win_args
{
	Call_winManager * call_winManager_p;
	vector<Soap_format> * alignment_vec_p;
	gzoutstream * consensus_p;
	my_ofstream * baseinfo_p;
	Genome * genome;
	Prob_matrix * mat;
	Parameter *para;
	int index;
	int ali_index;
	SfsMethod * sfsMethod_p;

	// initialize
	inline _call_win_args(Call_winManager* a, vector<Soap_format> * b,gzoutstream * c,my_ofstream* d,Genome* e,Prob_matrix* f,Parameter* g,int h, int i, SfsMethod *j)
	{
		call_winManager_p = a;
		alignment_vec_p = b;
		consensus_p = c;
		baseinfo_p = d;
		genome = e;
		mat = f;
		para = g;
		index = h;
		ali_index = i;
		sfsMethod_p = j;
	};
	inline _call_win_args()
	{
		call_winManager_p = NULL;
		alignment_vec_p = NULL;
		consensus_p = NULL;
		baseinfo_p = NULL;
		genome = NULL;
		mat = NULL;
		para = NULL;
		index = 0;
		ali_index = 0;
		sfsMethod_p = NULL;
	};
}CALL_WIN_ARGS;

//void* _call_winManager_soap2cns_norm(void * __args);

class Call_win_Task : public CTask
{
public:
	// overload the function Run().
	int Run();
	/*{
		_call_winManager_soap2cns_norm(this->m_ptrData);
		return 1;
	}*/
};
