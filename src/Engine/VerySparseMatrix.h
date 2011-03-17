// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
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

/*! \file VerySparseMatrix.h
 *
 *  A class to represent a sparse matrix in trivial format
 *
 */
 
#ifndef VERY_SPARSE_MATRIX_HEADER_H
#define VERY_SPARSE_MATRIX_HEADER_H


namespace Dmrg {
	//! Yet another sparse matrix class
	template<class T>
	class VerySparseMatrix {
		typedef std::pair<size_t,size_t> PairType;
	public:
		VerySparseMatrix(size_t rank) : rank_(rank),sorted_(true)
		{}

		template<typename CrsMatrixType>
		VerySparseMatrix(const CrsMatrixType& crs)
		{
			rank_=crs.rank();
			sorted_=true;
			for (size_t i=0;i<rank_;i++) {
				for (int k=crs.getRowPtr(i);k<crs.getRowPtr(i+1);k++) { 
					// (i,crs.getCol(k)) --> coordinate
					PairType coordinate = PairType(i,crs.getCol(k));
					coordinates_.push_back(coordinate);
					values_.push_back(crs.getValue(k));
				}
			}
		}

		VerySparseMatrix(const VerySparseMatrix& vsm,const T& eps)
		{
			rank_=vsm.rank();
			for (size_t i=0;i<vsm.values_.size();i++) {
				if (fabs(vsm.values_[i])<eps) continue;
				coordinates_.push_back(vsm.coordinates_[i]);
				values_.push_back(vsm.values_[i]);
			}
			sorted_=true;
		}

		T& operator()(size_t row,size_t col)
		{
			std::pair<size_t,size_t> coordinate(row,col);
			int x=utils::isInVector(coordinates_,coordinate);
			if (x<0) {
				coordinates_.push_back(coordinate);
				values_.push_back(0);
				x=values_.size()-1;
				sorted_=false;
			}
			return values_[x];
		}

		bool operator==(const VerySparseMatrix<T>& other) const 
		{
			if (rank_!=other.rank_) return notEqual("rank");

			if (!utils::vectorEqual(values_,other.values_)) return notEqual("values");
			if (!utils::vectorEqual(coordinates_,other.coordinates_)) return notEqual("coordinates");
			if (sorted_!=other.sorted_) return notEqual("sorted");
			return true;
		}

		bool operator!=(const VerySparseMatrix<T>& other) const 
		{
			return !(*this == other);
		}

		void clear()
		{
			rank_=0;
			values_.clear();
			coordinates_.clear();
			sorted_=false;
		}

		template<typename CrsMatrixType>
		void operator=(const CrsMatrixType& crs)
		{
			clear();
			rank=crs.rank();
			sorted_=true;
			for (size_t i=0;i<rank_;i++) 
				for (int k=crs.getRowPtr(i);k<crs.getRowPtr(i+1);k++) 
					// (i,crs.getCol(k)) --> coordinate
					this->operator()(i,crs.getCol(k))=crs.getValue(k);
					// crs.getValue(k) --> value			
		}

		//! same as T& operator() but doesn't check for dupes
		//! faster than T& operator() but use with care.
		void set(size_t row,size_t col,const T& value)
		{
			std::pair<size_t,size_t> coordinate(row,col);
			coordinates_.push_back(coordinate);
			values_.push_back(value);
		}

		T operator()(size_t row,size_t col) const
		{
			std::pair<size_t,size_t> coordinate(row,col);
			int x=utils::isInVector(coordinates_,coordinate);
			if (x<0) {
				//T value=0;
				return 0;
			}
			return values_[x];
		}

		const T& operator()(size_t i) const
		{
			return values_[i];
		}

		T& operator()(size_t i) 
		{
			return values_[i];
		}

		void operator+=(VerySparseMatrix<T>& other)
		{
			if (sorted_ && other.sorted()) plusEqualOrd(other);
			else throw std::runtime_error("VerySparseMatrix::operator+=(): unsorted\n"); //plusEqualUnord(other);
		}

		void sort()
		{
			std::vector<size_t> perm;
			sort(perm);
		}

		void sort(std::vector<size_t>& perm)
		{
			std::vector<size_t> rows(coordinates_.size()),cols(coordinates_.size());
			std::vector<T> values = values_;
			for (size_t i=0;i<coordinates_.size();i++) {
				rows[i]=coordinates_[i].first;
				cols[i]=coordinates_[i].second;
				
			}

			perm.resize(rows.size());
			utils::sort(rows,perm);

			for (size_t i=0;i<coordinates_.size();i++) {
				coordinates_[i].first = rows[i];
				coordinates_[i].second = cols[perm[i]];
				values_[i] = values[perm[i]];
			}
			sorted_=true;
		}

		size_t rank() const { return rank_; }

		size_t getRow(size_t i) const { return coordinates_[i].first; }

		void getRow(std::vector<size_t>& cols,size_t row,size_t startIndex=0) const
		{
			cols.clear();
			for (size_t i=startIndex;i<coordinates_.size();i++) {
				if (coordinates_[i].first == row) cols.push_back(i);
				if (sorted_ && coordinates_[i].first > row) break;
			}
		}

		size_t getColumn(size_t i) const { return coordinates_[i].second; }

		void getColumn(std::vector<size_t>& rows,size_t col) const
		{
			rows.clear();
			for (size_t i=0;i<coordinates_.size();i++) 
				if (coordinates_[i].second == col) rows.push_back(i);
		}

		size_t nonZero() const { return values_.size(); }

		T getValue(size_t i) const { return values_[i]; }

		template<typename T1>
		friend std::ostream& operator<<(std::ostream& os,const VerySparseMatrix<T1>& m);

		template<typename IoType>
		void saveToDisk(IoType& outHandle)
		{
			std::string s = "rank="+utils::ttos(rank_);
			outHandle.printline(s);
			outHandle.printVector(coordinates_,"coordinates");
			outHandle.printVector(values_,"values");
			s="#######\n";
			outHandle.printline(s);
		}

		template<typename IoType>
		void loadFromDisk(IoType& inHandle,bool check=false)
		{
			clear();
			std::string s = "rank=";
			inHandle.readline(rank_,s);
			inHandle.read(coordinates_,"coordinates");
			if (check) checkCoordinates();
			inHandle.read(values_,"values");
			sorted_=true;
		}

		//! for debuggin only:
		void checkCoordinates() const
		{
			size_t flag=0;
			for (size_t i=0;i<coordinates_.size();i++) {
				if (coordinates_[i].first<0 || 	coordinates_[i].first>=rank_) {
					std::cerr<<"coordinates["<<i<<"].first="<< coordinates_[i].first<<"\n";
 					flag=1;
					break;
				}
				if (coordinates_[i].second<0 || 	coordinates_[i].second>=rank_) {
					std::cerr<<"coordinates["<<i<<"].second="<< coordinates_[i].second<<"\n";
					flag=2;
					break;
				}
			}
			if (flag==0) return;
			std::cerr<<"rank="<<rank_<<"\n";
			throw std::runtime_error("VerySparseMatrix::checkCoordinates()\n");
		}
		
		bool sorted() const { return sorted_; }
		
		
	private:	
		size_t rank_;
		std::vector<T> values_;
		std::vector<PairType> coordinates_;
		bool sorted_;
		
		void plusEqualOrd(VerySparseMatrix<T>& other)
		{
			size_t j = 0;
			size_t i = 0;
			// pre-alloc memory:
			std::vector<PairType> newcoord(coordinates_.size()+other.coordinates_.size());
			std::vector<T> newvals(coordinates_.size()+other.coordinates_.size());
			size_t counter = 0;
			while(i<coordinates_.size() && j<other.coordinates_.size()) {
				PairType row = coordinates_[i];
				PairType row2 = other.coordinates_[j];
				//std::cerr<<row<<" "<<row2<<"\n";
				if (row2 < row) {
					//newcoord.push_back(j);
					newcoord[counter]=row2;
					newvals[counter]=other.values_[j];
					counter++;
					j++;
				} else if (row2 == row) {
					newcoord[counter]=row2;
					newvals[counter]=values_[i] + other.values_[j];
					j++;
					i++;
					counter++;
				} else {
					newcoord[counter]=row;
					newvals[counter]=values_[i];
					i++;
					counter++;
				}
			}
			for (size_t j1=j;j1<other.coordinates_.size();j1++) {
				newcoord[counter]=other.coordinates_[j1];
				newvals[counter]=other.values_[j1];
				counter++;
			}
			for (size_t j1=i;j1<coordinates_.size();j1++) {
				newcoord[counter]=coordinates_[j1];
				newvals[counter]=values_[j1];
				counter++;
			}
			coordinates_.resize(counter);
			values_.resize(counter);
			for (size_t j1=0;j1<counter;j1++) {
				coordinates_[j1] = newcoord[j1];
				values_[j1] = newvals[j1];
			}
			// this is sorted
		}
		
// 		void plusEqualUnord(VerySparseMatrix<T>& other)
// 		{
// 			throw std::runtime_error("dfjdkfdf\n");
// 			std::vector<size_t> cols1,cols2;
// 			size_t startIndex1=0,startIndex2=0;
// 			std::vector<PairType> newCoordinates;
// 			std::vector<T> newValues;
// 			
// 			for (size_t i=0;i<other.rank_;i++) {
// 				other.getRow(cols1,i,startIndex1);
// 				getRow(cols2,i,startIndex2);
// 				sumRows(newCoordinates,newValues,i,other,cols1,cols2);
// 			}
// 			for (size_t i=0;i<newValues.size();i++) {
// 				coordinates_.push_back(newCoordinates[i]);
// 				values_.push_back(newValues[i]);
// 			}
// 			if (coordinates_.size()>0) sort();
// 		}
		
// 		void sumRows(
// 			std::vector<PairType>& newCoordinates,
// 			std::vector<T>& newValues,
// 			size_t thisRow,
//    			const VerySparseMatrix<T>& other,
//       			const std::vector<size_t>& cols1,
// 	 		const std::vector<size_t>& cols2)
// 		{
// 			std::vector<size_t> realCols2(cols2.size());
// 			for (size_t i=0;i<cols2.size();i++) realCols2[i] = getColumn(cols2[i]);
// 
// 			for (size_t i=0;i<cols1.size();i++) {
// 				size_t col = other.getColumn(cols1[i]);
// 				int x = utils::isInVector(realCols2,col);
// 				if (x<0) {
// 					// add entry (thisRow,col) = other.getValue(cols1[i])
// 					newCoordinates.push_back(PairType(thisRow,col));
// 					newValues.push_back(other.getValue(cols1[i]));
// 					continue;
// 				}
// 				values_[cols2[x]] += other.getValue(cols1[i]);
// 			}
// 		}

		bool notEqual(const char *s) const
		{
			std::cerr<<"notEqual="<<s<<"\n";
			return false;
		}
	}; // VerySparseMatrix

	template<typename T>
	std::ostream& operator<<(std::ostream& os,const VerySparseMatrix<T>& m)
	{
		for (size_t i=0;i<m.values_.size();i++) 
			os<<"verysparse("<<m.coordinates_[i].first<<","<<m.coordinates_[i].second<<")="<<m.values_[i]<<"\n";
		return os;
	}

	template<typename T>
	bool isHermitian(const VerySparseMatrix<T>& m)
	{
		T eps =1e-6;
		for (size_t i=0;i<m.nonZero();i++) { 
			size_t row = m.getRow(i);
			size_t col = m.getColumn(i);
			if (fabs(m.getValue(i)-m(col,row))>eps) {
				return false;
			}
		}
		return true;
	}
} // namespace Dmrg
/*@}*/	
#endif
