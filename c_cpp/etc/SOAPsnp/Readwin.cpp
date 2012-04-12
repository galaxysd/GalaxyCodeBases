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

#include "Readwin.h"


Readwin::Readwin(void)
	: m_win_idx(0)
	, m_win_pos(0)
	, m_win_size(0)
	, m_last_pos(0)
{
	m_ail_vec[0].clear();
	m_ail_vec[1].clear();
	m_soap2cns_idx[0] = 0;
	m_soap2cns_idx[1] = 0;
	m_win_size = WIN_SIZE;
	m_win_pos -= m_win_size;
}

int Readwin::m_last_count = 0;

Readwin::~Readwin(void)
{
}

/**
 * DATE: 2010-8-9
 * FUNCTION: return the m_soap2cns_idx
 * PARAMETER: read: a record which will be add to the Soap_format vector.
 * RETURN:	READ_POS_ERROR when read's position less than first win's position.
 *			READ_EXCEED_WIN when read's position exceed the first win's tail.
 *			ADD_SUCCESSFUL when read add successful.
 *			COME_NEW_CHR when the record is coming into new chromosome.
 */
int Readwin::addRead(const std::string& read)
{
	std::istringstream ss(read);
	Soap_format soap;
	if (ss >> soap)
	{
		if (soap.get_pos() < m_win_pos)
		{
			// if the new record is not the same with the last record, then set a new win.
			if (soap.get_chr_name() != m_chr_name)
			{
				m_last_pos = m_win_pos + 2 * m_win_size;
				m_last_count++;	// add up the number of the Readwin that reach pre-chromosome's last.
				m_last_record = soap;	// keep the last record.
				return COME_NEW_CHR;
			}
			// the read's position was not sorted when it's position less than first win's position.
			return READ_POS_ERROR;
		}

		// add the soap object to the first win when it's pos between the first win.
		m_chr_name = soap.get_chr_name();

		if (soap.get_pos() < (m_win_pos + m_win_size))
		{
			m_ail_vec[m_win_idx].push_back(soap);
			if ((soap.get_pos() + soap.get_read_len()) >= (m_win_pos + m_win_size))
			{
				m_ail_vec[1 - m_win_idx].push_back(soap);
			}
		}
		else if (soap.get_pos() < (m_win_pos + 2 * m_win_size))
		{
			// add the soap object to the second win when it's pos between the second win.
			m_ail_vec[1 - m_win_idx].push_back(soap);
			// set begin index of second win
			m_soap2cns_idx[1 - m_win_idx] = m_ail_vec[1 - m_win_idx].size() - 1;
			return READ_EXCEED_WIN;
		}
		else
		{
			// store the soap object when it's position exceed the second win.
			m_last_record = soap;
			// set begin index of second win
			if (m_ail_vec[1 - m_win_idx].size() != 0)
			{
				m_soap2cns_idx[1 - m_win_idx] = m_ail_vec[1 - m_win_idx].size();
			}
			else
			{
				m_soap2cns_idx[1 - m_win_idx] = 0;
			}
			return READ_EXCEED_WIN;
		}
		return ADD_SUCCESSFUL;
	}
	return RECORD_ERROR;
}

/**
 * DATE: 2010-8-9
 * FUNCTION: return the m_soap2cns_idx
 * PARAMETER: void
 * RETURN: the index of m_ali_vec that can be processed in the beginning of the 'soap2cns'
 */
int Readwin::getSoap2cnsIdx(void)
{
	return m_soap2cns_idx[m_win_idx];
}

/**
 * DATE: 2010-8-9
 * FUNCTION: return the read win that can be processed
 * PARAMETER: void
 * RETURN: the win that can be processed.
 */
vector<Soap_format> &Readwin::getReadwin(void)
{
	return m_ail_vec[m_win_idx];
}

/**
 * DATE: 2010-8-9
 * FUNCTION: judge if m_last_record is empty.
 * PARAMETER: pos: the last record's pos.
 * RETURN: if m_last_record is empty.
 */
bool Readwin::isAbleToAdd()
{

	if (m_last_pos != 0)
	{
		// because the last record reach the last read of pre-chromosome, so cann't be added.
		return false;
	}


	if (m_last_record.get_pos() < (m_win_pos + m_win_size))
	{
		if ((m_last_record.get_pos() + m_last_record.get_read_len()) >= (m_win_pos + m_win_size))
		{
			// add it to the second win.
			m_ail_vec[1 - m_win_idx].push_back(m_last_record);
		}
		return true;
	}
	else if (m_last_record.get_pos() < (m_win_pos + 2 * m_win_size))
	{
		// if last record inside second win, then add it to the second win.
		m_ail_vec[1 - m_win_idx].push_back(m_last_record);
		m_soap2cns_idx[1 - m_win_idx] = m_ail_vec[1 - m_win_idx].size() - 1;
	}
	return false;
}


/**
 * DATE: 2010-8-9
 * FUNCTION: release the win that was processed, and move it behind anorther win.
 * PARAMETER: void
 * RETURN: void.
 */
void Readwin::winChange(void)
{
	// clear the vector had been processed.
	m_ail_vec[m_win_idx].clear();
	// move behind another win.
	m_win_idx = 1 - m_win_idx;
	// position move to another win.
	m_win_pos += m_win_size;
}

/**
 * DATE: 2010-8-26
 * FUNCTION: return the m_last_count
 * PARAMETER: void
 * RETURN: m_last_count.
 */
const int Readwin::getLastCount(void)
{
	return m_last_count;
}

/**
 * DATE: 2010-8-26
 * FUNCTION: reset the Readwin's data.
 * PARAMETER: void
 * RETURN: 0
 */
int Readwin::reset(void)
{
	m_win_pos = 0 - m_win_size;
	m_last_pos = 0;
	m_last_count = 0;
	m_soap2cns_idx[0] = 0;
	m_soap2cns_idx[1] = 0;
	return 0;
}


// set the m_last_pos's value.
/**
 * DATE: 2010-9-3
 * FUNCTION: set the m_last_pos's value.
 * PARAMETER: void
 * RETURN: void
 */
void Readwin::setLastPos(int last_pos)
{
	m_last_pos = last_pos;
}
