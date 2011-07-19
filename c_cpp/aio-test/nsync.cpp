/*
 * nsync.cpp
 *
 *  Created on: Jul 5, 2010
 *      Author: Aqua
 */

#include"nsync.h"

using namespace std;
using namespace boost;

namespace bp = ::boost::process;
namespace ba = ::boost::asio;

#define BUFFER_SIZE 0xFFFFFF

int ImportFileList(nSyncBuf_s & buf, char * fn)
{
	char errStr[4096];
	int count(0);

	if(access(fn, R_OK) == -1)
	{
		snprintf(errStr, 4095, "[%s][%s][%d]: File->%s", __FILE__, __FUNCTION__, __LINE__, fn);
		perror(errStr);
		return -1;
	}
	ifstream I(fn);
	string strTmp;
	while(true)
	{
		getline(I, strTmp, '\n');
		if(!I)
			break;
		if(access(strTmp.c_str(), R_OK) == -1)
		{
			snprintf(errStr, 4095, "[%s][%s][%d]: File->%s", __FILE__, __FUNCTION__, __LINE__, strTmp.c_str());
			perror(errStr);
			continue;
		}
		buf.fnList.push_back(strTmp);
		++count;
	}

	return count;
}

int InsertFile(nSyncBuf_s & buf, char * fn)
{
	char errStr[4096];

	if(access(fn, R_OK) == -1)
	{
		snprintf(errStr, 4095, "[%s][%s][%d]: File->%s", __FILE__, __FUNCTION__, __LINE__, fn);
		perror(errStr);
		return -1;
	}
	else
		buf.fnList.push_back(fn);
	return 1;
}

bp::child StartChild(nSyncBuf_s & buf)
{
	string exec = bp::find_executable_in_path("zcat");

	vector<string> args;
	args.push_back(bp::find_executable_in_path("zcat"));
	args.push_back("-f");

	for(uint i(0); i < buf.fnList.size(); ++i)
	{
		args.push_back(buf.fnList[i]);
	}

	bp::context ctx;
	ctx.stdout_behavior = bp::capture_stream();
	ctx.stderr_behavior = bp::inherit_stream();
	ctx.environment = bp::self::get_environment();

	return bp::launch(exec, args, ctx);
}

void endRead(nSyncBuf_s &, const boost::system::error_code & ec, size_t bytes_transferred);

void nSyncRead(nSyncBuf_s & buf)
{
	buf.in.async_read_some(boost::asio::buffer(buf.buffer, BUFFER_SIZE),
		boost::bind(&endRead, boost::ref(buf), ba::placeholders::error, ba::placeholders::bytes_transferred));
}

void endRead(nSyncBuf_s & buf, const boost::system::error_code & ec, size_t bytes_transferred)
{
	if(!ec)
	{
		if(buf.atomEnabler)
			buf.lock();
		buf.Routine(buf, bytes_transferred);
		nSyncRead(buf);
	}
}

void nSync(nSyncBuf_s & buf)
{
	if(buf.fnList.empty())
	{
		perror("Empty filename vector. Program Terminated");
		exit(EXIT_FAILURE);
	}
	bp::child c = StartChild(buf);

	buf.in.assign(c.get_stdout().handle().release());

	nSyncRead(buf);

	buf.io_service.run();

	c.wait();
}

nSyncBuf_s::nSyncBuf_s(_ConsumeBufferRoutine_f func, bool atomic) : atomEnabler(atomic), Routine(func), in(io_service)
{
	buffer = new char[BUFFER_SIZE];
	if(!buffer)
	{
		perror("[nSyncBuf_s]: Buffer allocation error.");
		throw bad_alloc();
	}
}

nSyncBuf_s::~nSyncBuf_s()
{
	if(buffer)
	{
		delete[] buffer;
	}
	else
	{
		perror("[nSyncBuf_s]: Null ptr of buffer when destruction.");
		throw bad_alloc();
	}
}

inline void nSyncBuf_s::lock()
{
	atom.lock();
}

inline void nSyncBuf_s::unlock()
{
	atom.unlock();
}

void nSyncBuf_s::release()
{
	atom.unlock();
}

nSyncBufAry_t nSyncBuf_s::ptr()
{
	return buffer;
}

nSyncBufAry_t nSyncBuf_s::operator *()
{
	return buffer;
}

