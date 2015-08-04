/*
 ******************************************************************************
 *Copyright 2010
 *	BGI-SHENZHEN
 *All Rights Reserved
 *ATHUOR : Bill Tang
 *CREATE DATE : 2010-8-9
 *CLASS NAME: Readwin
 *FUNCTION : 
 *FILE NAME : Readwin.cpp
 *UPDATE DATE : 2010-9-3
 *UPDATE BY : Bill Tang
 *
 *UPDATE: 2010-9-3 add function setLastPos.
 *******************************************************************************
 */

#pragma once
#include "soap_snp.h"
#include <pthread.h>

#define WIN_SIZE global_win_size // the win's sizes
#define READ_POS_ERROR -1	// the reads were not sorted.
#define ADD_SUCCESSFUL 1	// add read successful.
#define READ_EXCEED_WIN 2	// read's position exceed first win's position.
#define RECORD_ERROR -2		// the read record is error.
#define COME_NEW_CHR 3		// when the record is coming into new chromosome

#define BE_PROCESSED	-1	// the read win had been processed.

class Readwin
{
public:
	Readwin(void);
	virtual ~Readwin(void);
	virtual int addRead(const std::string& read);
private:
	// the index of m_ali_vec that can be processed in the beginning of the 'soap2cns'
	int m_soap2cns_idx[2];
	// the index of m_ali_vec[] to be processed.
	int m_win_idx;
	// the previous win's position
	int m_win_pos;
	// the win's size
	int m_win_size;
	// the two win vector
	vector<Soap_format> m_ail_vec[3];
	// the last record which is contained in neither m_ail_vec[0] nor m_ail_vec[1]
	Soap_format m_last_record;
	// record the position once one of the samples reach a break point that the chromosome changed.
	int m_last_pos;
	// record the number of wins which reached a chromosome's last position.
	static int m_last_count;
	// the thread synchronization lock.
	static pthread_mutex_t m_pthreadMutex;
public:
	// return the m_soap2cns_idx
	virtual int getSoap2cnsIdx(void);
	// return the read win that can be processed
	virtual vector<Soap_format> &getReadwin(void);
	// judge if m_last_record is empty.
	virtual bool isAbleToAdd();
	// release the win that was processed, and move it behind anorther win
	virtual void winChange(void);
	// return the m_last_count
	virtual int getLastCount(void);
	// reset the Readwin's data
	virtual int reset(void);
private:
	// store the current processed chromosome name
	std::string m_chr_name;
public:
	// set the m_last_pos's value.
	virtual void setLastPos(const int last_pos);
private:
	// record the m_ail_vec's index which to be processed.
	int m_process_idx;
	// wait for the read win to be proessed.
	virtual void setProcessIdx(void);
public:
	// set the start window position
	void setWin_pos(int start_pos);
};

