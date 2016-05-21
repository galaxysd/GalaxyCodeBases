#include<vector>
#include<unistd.h>
#include<pthread.h>
#include<sys/time.h>

using namespace std;

#ifndef _THREADMANAGER_H_
#define _THREADMANAGER_H_

typedef void* (*pthreadproc)(void*);

class CThreadManager
{
public:
	CThreadManager();
	CThreadManager(int cpu);
	~CThreadManager();

	void				AddThread(pthreadproc pProc, void* param);
	void				Run();
	void				SetTimer(float t = 0.5);
	void				Reset();
	void				SetCPU(int count);
private:
public:
	vector<pthreadproc>		m_vproc;
	vector<void*>			m_vparam;
	int				m_cpu;
	int				m_nmaxthread;
	int				m_nindex;

protected:
private:
	float				m_ftimer;
};

// ThreadParam
struct THRPARAM
{
	CThreadManager*		pMgr;
	int					flag;
};

#endif

