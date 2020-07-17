#ifndef FILE_REAL_EXCEPTIONS_H
#define FILE_REAL_EXCEPTIONS_H

namespace RealLib {

class RealLibException : public std::exception {
   char m_what[128];
public:
   RealLibException(const char *what = NULL) throw();
   virtual const char *what() const throw()
      {  return m_what; }
   virtual const char *kind() const throw()
	  {  return "RealLibException"; }
};

class PrecisionException : public RealLibException {
public:
   PrecisionException(const char *what = NULL) throw()
      : RealLibException(what) {}
   const char *kind() const throw()
	  {  return "PrecisionException"; }
};

class DomainException : public RealLibException {
public:
   DomainException(const char *what = NULL) throw()
      : RealLibException(what) {}
   virtual const char *kind() const throw()
	  {  return "DomainException"; }
};

static inline
std::ostream & operator << (std::ostream &os, RealLibException &e)
{
	return os << e.kind() << ": " << e.what();
}

}

#endif // FILE
