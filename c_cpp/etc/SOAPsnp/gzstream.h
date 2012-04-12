// ============================================================================
// gzstream, C++ iostream classes wrapping the zlib compression library.
// Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// File          : gzstream.h
// Revision      : $Revision: 1.5 $
// Revision_date : $Date: 2002/04/26 23:30:15 $
// Author(s)     : Deepak Bandyopadhyay, Lutz Kettner
// 
// Standard streambuf implementation following Nicolai Josuttis, "The 
// Standard C++ Library".
// ============================================================================

#ifndef GZSTREAM_H
#define GZSTREAM_H 1

// standard C++ with new header file names and std:: namespace
#include <iostream>
#include <fstream>
//#include <ostream>
#include <zlib.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#ifdef GZSTREAM_NAMESPACE
namespace GZSTREAM_NAMESPACE {
#endif

// ----------------------------------------------------------------------------
// Internal classes to implement gzstream. See below for user classes.
// ----------------------------------------------------------------------------

class gzstreambuf : public std::streambuf {
private:
    static const int bufferSize = 47+256;    // size of data buff
    // totals 512 bytes under g++ for igzstream at the end.

    gzFile           file;               // file handle for compressed file
    char             buffer[bufferSize]; // data buffer
    int             opened;             // open/close state of stream
    int              mode;               // I/O mode

    int flush_buffer();
public:
    gzstreambuf() : opened(0) 
	{
        setp( buffer, buffer + (bufferSize-1));
        setg( buffer + 4,     // beginning of putback area
              buffer + 4,     // read position
              buffer + 4);    // end position      
        // ASSERT: both input & output capabilities will not be used together
    }
   int is_open() 
	{ 
		return opened;
	}
    gzstreambuf* open( const char* name, int open_mode);
    gzstreambuf* close();
    ~gzstreambuf() { close(); }
    
    virtual int     overflow( int c = EOF);
    virtual int     underflow();
    virtual int     sync();
	gzstreambuf& operator=(const gzstreambuf& gzs);
};

class gzstreambase : virtual public std::ios 
{
protected:
    gzstreambuf buf;
public:
    gzstreambase() 
	{ 
		init(&buf);
	}
    gzstreambase( const char* name, int open_mode);
    ~gzstreambase();
	int is_open() 
	{ 
		return buf.is_open(); 
	}
    void open( const char* name, int open_mode);
    void close();
    gzstreambuf* rdbuf() 
	{ 
		return &buf;
	}
	gzstreambase& operator=(const gzstreambase& gzst);
};

// ----------------------------------------------------------------------------
// User classes. Use igzstream and ogzstream analogously to ifstream and
// ofstream respectively. They read and write files based on the gz* 
// function interface of the zlib. Files are compatible with gzip compression.
// ----------------------------------------------------------------------------

class igzstream : public gzstreambase, public std::istream 
{
public:
    igzstream() : std::istream( &buf) 
	{
		/*gzstreambase::init(&buf);*/
	} 
    igzstream( const char* name, int open_mode = std::ios::in)
        : gzstreambase( name, open_mode), std::istream( &buf) 
	{
		/*gzstreambase::init(&buf);*/
	}  
    gzstreambuf* rdbuf() { return gzstreambase::rdbuf(); }
    void open( const char* name, int open_mode = std::ios::in)
	{
        gzstreambase::open( name, open_mode);
    }
};

class ogzstream : public gzstreambase, public std::ostream
{
public:
    ogzstream() : std::ostream( &buf) 
	{
		/*gzstreambase::init(&buf);*/
	}
    ogzstream( const char* name, int mode = std::ios::out)
        : gzstreambase( name, mode), std::ostream( &buf) 
	{
		/*gzstreambase::init(&buf);*/
	}  
    gzstreambuf* rdbuf() 
	{ 
		return gzstreambase::rdbuf(); 
	}
    void open( const char* name, int open_mode = std::ios::out)
	{
        gzstreambase::open( name, open_mode);
    }
	ogzstream& operator=(const ogzstream& ogzst);
};


class myfstream
{
	public:
		virtual void open(const char* name, int open_mode = std::ios::out) = 0;
		virtual void close() = 0;
		virtual void clear() = 0;
		virtual bool good() =0;
		virtual std::ostream& flush() =0;
		virtual std::ostream& endl() =0;
		virtual void write(const char* name, int streamsizecount) =0;
		virtual int is_open() = 0;
		virtual myfstream & operator<<(const double value) = 0;
		virtual myfstream & operator<<(const char* value) = 0;
		virtual myfstream & operator<<(const char value) = 0;
		virtual myfstream & operator<<(const unsigned char value) = 0;
		virtual myfstream & operator<<(const std::string & value) = 0;
		virtual myfstream & operator<<(const int value) = 0;
		virtual myfstream & operator<<(const int* value) = 0;
		virtual myfstream & operator<<(const unsigned int value) = 0;
		virtual myfstream & operator<<(const float value) = 0;
		virtual int set_f(std::ios_base::fmtflags value) = 0;
		//virtual ostream & operator<<(std::ostream& value) = 0;
		//virtual myfstream & operator<< (std::ostream& (__cdecl *fun) (ostream& value)) =0;

};

/**
 * a class to use ofstream.
 *
 */
class myofstream : public myfstream
{
	private:
		std::ofstream ogz;
	public:
		myofstream() : ogz(){}
		myofstream(const char * name, int mode = std::ios::out)
		       :ogz(name){}
		virtual inline void open( const char* name, int open_mode = std::ios::out) 
		{
			ogz.open(name);
		}
		virtual inline void close()
		{
			ogz.close();
		}
		virtual inline void write(const char* name, int streamsizecount )
		{			
			ogz.write(name, streamsizecount);
		}
		virtual inline std::ostream& flush()
		{
			return ogz.flush();
		}
		virtual inline std::ostream& endl()
		{
			ogz << "\n";
			return ogz.flush();
		}
		virtual inline void clear()
		{
			ogz.clear();
		}
		virtual inline bool good()
		{
			return ogz.good();
		}
		virtual inline int is_open()
		{
			return ogz.is_open();
		}
		
		virtual inline myofstream & operator<<(const double value)
		{
			ogz << value;			
			return *this;
		}
		virtual inline myofstream & operator<<(const char * value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myofstream & operator<<(const int value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myofstream & operator<<(const int* value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myofstream & operator<<(const unsigned int value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myofstream & operator<<(const char value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myofstream & operator<<(const unsigned char value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myofstream & operator<<(const std::string & value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myofstream & operator<<(const float value)
		{
			ogz << value;
			return *this;
		}
		virtual inline int set_f(std::ios_base::fmtflags value)
		{
			return ogz.setf(value);
		}
		/*virtual inline myofstream & operator<<(const std::ios& value)
		{
			ogz << value;
			return *this;
		}*/
		/*virtual myofstream& operator<< (std::ostream& (__cdecl *fun)(ostream& value))
		{
			(*fun)(*this);
			return *this;
		}*/

};

/**
 * a class to use ogzstream.
 *
 */
class myogzstream : public myfstream
{
	private:
		ogzstream ogz;
	public:
		myogzstream() : ogz(){}
		myogzstream(const char * name, int mode = std::ios::out)
		       :ogz(name, mode){}
		virtual inline void open( const char* name, int open_mode = std::ios::out) 
		{
			ogz.open(name, open_mode);			
		}
		virtual inline void close()
		{
			ogz.close();
		}
		virtual inline bool good()
		{
			return ogz.good();
		}
		virtual inline void write(const char* name, int streamsizecount)
		{			
			ogz.write(name,streamsizecount);
		}
		virtual inline void clear()
		{
			ogz.clear();
		}
		virtual std::ostream& flush()
		{
			return ogz.flush();
		}
		virtual std::ostream& endl()
		{
			ogz<<"\n";
			return ogz.flush();
		}
		virtual inline int is_open()
		{
			return ogz.is_open();
		}
		
		virtual inline myogzstream & operator<<(const double value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myogzstream & operator<<(const char * value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myogzstream & operator<<(const int value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myogzstream & operator<<(const int* value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myogzstream & operator<<(const unsigned int value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myogzstream & operator<<(const char value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myogzstream & operator<<(const unsigned char value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myogzstream & operator<<(const std::string & value)
		{
			ogz << value;
			return *this;
		}
		virtual inline myogzstream & operator<<(const float value)
		{
			ogz << value;
			return *this;
		}
		virtual inline int set_f(std::ios_base::fmtflags value)
		{
			return ogz.setf(value);
		}
		/*virtual inline myogzstream & operator<<(const std::ios& value)
		{
			ogz << value;
			return *this;
		}*/
		/*virtual myogzstream& operator<< (std::ostream& (__cdecl *fun)(ostream& value))
		{
			(*fun)(*this);
			return *this;
		}*/
};

#ifdef GZSTREAM_NAMESPACE
} // namespace GZSTREAM_NAMESPACE
#endif

#endif // GZSTREAM_H
// ============================================================================
// EOF //


