/*
 * Asynchronous IO buffer adapter
 *
 * This library is under no license, hence,
 * you can do whatever you wanna do on this
 * library.
 *
 *  Created on: Jul 5, 2010
 *      Author: Aqua
 */

#ifndef _NSYNC_H_AQUA
#define _NSYNC_H_AQUA

/*
 * Here is a illustration on you use this lib,
 * this lib requires you to define a Routine that
 * flush the data from buffer and do what you wanna
 * ----------------------------------------------------------------------
	void Routine(nSyncBuf_s & buf, size_t size)
	{
		cout<<string(*buf, size);
		buf.release();
	}

	int main(int argc, char ** argv)
	{
		nSyncBuf_s buf(Routine, true);
		cerr<<"[ImportFileList]:"<<ImportFileList(buf, argv[1])<<endl;
		nSync(buf);
	}
 * ----------------------------------------------------------------------
	nSyncBuf_s buf(coutRoutine, true):
		"nSyncBuf_s" is a class that contains the buffer
		and linked to boost::asio for asynchronous stream.
		"Routine" is the Routine that handles the data.
		You can pass the second parameter as "true" if the
		next run of "Routine" would depend on the accomplishment
		of the last run, thus to say, read and process sequentially,
		or "false" if "Routine"s are absolutely independent.
		PLEASE NOTE THE, for passing "true", it's your duty to
		"release" when the previous Routine was done. You can
		see "buf.release();" after "cout<<string(*buf, size);".
		Deadlock would occur if you forget to release when "true";
	cout<<string(*buf, size):
		This demonstrated the way to visit the buffer, "*buf"
		as well as "buf.ptr()" would return the foremost
		pointer of buffer. "size" indicates the length of
		effective data in buffer.
	ImportFileList(buf, argv[1]):
		"argv[1]" should be a filename that contains filenames,
		one row each filename.
	InsertFile(buf, argv[1]):
		"argv[1]" should be a single filename pending to read.
		Please notice that, for multiple files, the would be no
		EOF at the end of each file, thus, the data of next file
		(if available) would append tightly at the very end of
		previous file.
	nSync(buf):
		Let's rock.
 */


#if defined(__CYGWIN__)
#  define _WIN32_WINNT 0x0501
#  define __USE_W32_SOCKETS
#  undef BOOST_POSIX_API
#  define BOOST_WINDOWS_API
#endif
#include <boost/asio.hpp>
#define BOOST_PROCESS_WINDOWS_USE_NAMED_PIPE
#include <boost/process.hpp>
#include <boost/bind.hpp>
#include <unistd.h>
#include <fstream>
#include <boost/thread/mutex.hpp>
#include <exception>
#include <vector>

using namespace std;
using namespace boost;
namespace bp = ::boost::process;
namespace ba = ::boost::asio;

typedef char buf_t;
typedef buf_t * nSyncBufAry_t;

class nSyncBuf_s;
typedef void(*_ConsumeBufferRoutine_f)(nSyncBuf_s &, size_t);

class nSyncBuf_s
{
public:

	explicit nSyncBuf_s(_ConsumeBufferRoutine_f, bool);
	~nSyncBuf_s();

	void lock() __attribute__((always_inline));
	void unlock() __attribute__((always_inline));
	void release();
	nSyncBufAry_t ptr();
	nSyncBufAry_t operator *();

private:
	friend int ImportFileList(nSyncBuf_s &, char *);
	friend int InsertFile(nSyncBuf_s &, char *);
	friend void nSync(nSyncBuf_s &);
	friend bp::child StartChild(nSyncBuf_s &);
	friend void nSyncRead(nSyncBuf_s &);
	friend void endRead(nSyncBuf_s &, const boost::system::error_code &, size_t);

protected:
	vector<string> fnList;
	nSyncBufAry_t buffer;
	bool atomEnabler;
	mutex atom;
	_ConsumeBufferRoutine_f Routine;

public:
	//ASIO stream connectors
	ba::io_service io_service;
	#if defined(BOOST_POSIX_API)
	ba::posix::stream_descriptor in;
	#elif defined(BOOST_WINDOWS_API)
	ba::windows::stream_handle in;
	#else
	#error Unsupported platform.
	#endif
};

int ImportFileList(nSyncBuf_s &, char *) __attribute__((warn_unused_result));
int InsertFile(nSyncBuf_s &, char *) __attribute__((warn_unused_result));
void nSync(nSyncBuf_s &);

#endif /* _NSYNC_H_AQUA */
