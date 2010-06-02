// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."
 
*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************


*/
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file SimpleReader.h
 *
 *  Parse a unix file
 *
 */
#ifndef SIMPLE_READER_H
#define SIMPLE_READER_H
#include "Utils.h"


namespace Dmrg {
	class SimpleReader {
	public:
		SimpleReader(const std::string& fileName): fileName_(fileName),fin_(fileName.c_str())
		{
			if (!fin_ || !fin_.good() || fin_.bad()) throw std::runtime_error("Cannot open\n"); 
		}

		~SimpleReader() { fin_.close(); }

		template<typename ParametersType>
		void load(ParametersType& parameters) 
		{
			parameters <= *this;
		}

		void read(std::string& s)  { fin_>>label_; fin_>>s; }

		template<typename FieldType>
		void read(std::vector<FieldType>& v) 
		{
			// read label
			fin_>>label_;
			// read vector length
			size_t size;
			fin_>>size;
			// read vector
			v.resize(size);
			for (size_t i=0;i<v.size();i++) {
				fin_>>v[i];
			}
		}

		template<typename FieldType>
		void read(psimag::Matrix<FieldType>& v) 
		{
			// read label
			fin_>>label_;
			// read row and col
			size_t row,col;
			fin_>>row;
			fin_>>col;
			// read matrix
			v.resize(row,col);
			for (size_t i=0;i<v.n_row();i++) 
				for (size_t j=0;j<v.n_col();j++) 
					fin_>>v(i,j);
		}

		template<typename FieldType>
		void read(FieldType& x)  { fin_>>label_; fin_>>x; }

	private:
		std::string    fileName_,label_;
		std::ifstream fin_;
	};
} // namespace Dmrg

/*@}*/
#endif


