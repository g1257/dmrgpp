//-*-C++-*-

/** \ingroup MatrixParser */
/*@{*/

/*! \file matrixParser.h  
 *
 *  
 *
 */

#ifndef  JsonParser_MatrixParser_H
#define  JsonParser_MatrixParser_H

#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include "CharacterMapper.h"
#include "MatrixLike.h"

namespace JsonParser {
  
  template<typename MatrixLikeType>
  class MatrixParser: public CharacterMapper {
  public:

    typedef typename MatrixLikeType::value_type FieldType;

    MatrixLikeType& mat;
    std::wifstream& inputStream;
    int endPos;

    MatrixParser(std::wifstream& is, 
		 int end,
		 MatrixLikeType& m):
      mat(m),
      inputStream(is),
      endPos(end)
    {
      consumeMatrix();
    }
    
    //======================================================================

    bool consume(wchar_t c) {
      wchar_t nextChar(0);
      while(true) {
	nextChar = inputStream.get();
	if (nextChar == c)       return true;
	if (!inputStream.good()) return false;
	if (!isWhiteSpace(wcharToClass(nextChar)) )
	  break;
      }
      inputStream.unget();
      return false;
    }

    //======================================================================

    template<typename T>
    void consume(T& value) {
      inputStream >> value;
    }

    //======================================================================
    
    template<typename T>
    void consume(std::complex<T>& cmplx) {
      T r,i;
      if (!consume(L'['))
	throw std::logic_error("parsing error in consume(std::complex<Fieldtype> cmplx)\n  missing [\n");
      inputStream >> r;
      if (!consume(L','))
	throw std::logic_error("parsing error in consume(std::complex<Fieldtype> cmplx)\n  missing ,\n");
      inputStream >> i;
      if(!consume(L']'))
	throw std::logic_error("parsing error in consume(std::complex<Fieldtype> cmplx)\n  missing ]\n");
      cmplx = std::complex<T>(r,i);
    }
    
    //======================================================================
    
    void consumeMatrix()  {
      
      //Consume the lead bracket for the array
      if (!consume(L'['))
	throw std::logic_error("Missing [ marking the start of the array in consumeMatrix()\n");
	
      size_t datarow(0);

      while( !consume(L']')           && 
	     inputStream.good()       && 
	     inputStream.tellg() <= endPos) { // while we have not encountered the end of the array, etc.

	//Consume the lead bracket for the next row
	if(!consume(L'['))
	  throw std::logic_error("Missing [ marking the start of a row in consumeMarix()");
	
	size_t datacol(0);

	
	while(!consume(L']')) { // while we don't hit the end of the row 
	  
	  FieldType& value(mat(datarow,datacol));
	  consume(value);
          //std::cout << "mat(" << datarow << "," << datacol << ")= " << value << "\n";
	  datacol++;
	  consume(L',');

	}
	datarow++;
	consume(L',');
      }
      //      std::cout << "consume done\n";
      //      psimag::MatrixLike::printList(mat,std::cout);
    }
    
    //======================================================================

  };
  
} // end namespace JsonParser


/*@}*/
#endif

