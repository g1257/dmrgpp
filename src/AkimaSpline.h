#ifndef AKIMA_H_
#define AKIMA_H_

#include <cmath>
#include <stdexcept>
#include <utility>

//! A class to interpolate using akima spline
template<typename VectorType>
class AkimaSpline {
	typedef typename VectorType::ValueType RealType;
	
	struct AkimaStruct {
		RealType x0,x1,a0,a1,a2,a3;
	};
	
	
	public:
		typedef std::pair<RealType,RealType> IntervalType;

		AkimaSpline(const VectorType& x,const VectorType& s)
		{
			if (x.size()!=s.size())
				throw std::runtime_error("Number of X and Y points must be the same\n");
			if (x.size()<3) throw
				std::runtime_error("Number of X too small\n");
			VectorType sprime;
			calculateSprime(sprime,x,s);
			for (size_t i=0;i<x.size()-1;i++) {
				RealType u = x[i+1]-x[i];
				RealType a2 = u*(sprime[i+1]-sprime[i])-2*(s[i+1]-s[i]);
				RealType a3 = 3*u*(s[i+1]-s[i])-u*u*(sprime[i+1]+2*sprime[i]);
				RealType ucube = u*u*u;
				AkimaStruct ak;
				ak.x0=x[i]; ak.x1=x[i+1];
				ak.a0=s[i]; ak.a1=sprime[i]; ak.a2=a2/ucube; ak.a3=a3/ucube;
				akimaStruct_.push_back(ak);
			}
			size_t k = akimaStruct_.size()-1;
			interval_=IntervalType(akimaStruct_[0].x0,akimaStruct_[k].x1);
		}
			
		RealType operator()(const RealType& x) const
		{
			size_t i =findRange(x);
			return functionFor(x,i);
		}

		IntervalType getInterval() const { return interval_; }

	private:
		
		RealType functionFor(const RealType& x,size_t i) const
		{
			const AkimaStruct& ak = akimaStruct_[i];
			RealType u = x - ak.x0;
			return ak.a0+ak.a1*u+ak.a2*u*u+ak.a3*u*u*u;
		}
		
		// FIXME: make this function faster
		size_t findRange(const RealType& x) const
		{
			if (x<interval_.first || x>interval_.second)
				throw std::runtime_error("Akima: out of range\n");
			for (size_t i=0;i<akimaStruct_.size();i++)
				if (x>=akimaStruct_[i].x0 && x<=akimaStruct_[i].x1) return i;
	
			throw std::runtime_error("Akima: findRange(): internal error\n");
		}
		
		void calculateSprime(VectorType& sprime,const VectorType& x,const VectorType& s) const
		{
			VectorType d,w;
			calculateD(d,x,s);
			calculateW(w,d);
			sprime.resize(x.size());
			for (size_t i=0;i<sprime.size();i++) {
				RealType denom = w[i]+w[i+2];
				if (denom==0) {
					sprime[i] = (d[i+1]+d[i+2])*0.5;
					continue;
				}
				sprime[i] = w[i+2]*d[i+1]+w[i]*d[i+2];
				sprime[i] /= denom;
			}
		}
		
		void calculateD(VectorType& d,const VectorType& x,const VectorType& s) const
		{
			size_t k = x.size();
			d.resize(k+3);
			for (size_t i=2;i<=k;i++) 
				d[i] = (s[i-1]-s[i-2])/(x[i-1]-x[i-2]);
			d[1] = 2*d[2] - d[3];
			d[0] = 2*d[1] - d[2];
			d[k+1] = 2*d[k] - d[k-1];
			d[k+2] = 2*d[k+1]-d[k];
		}

		void calculateW(VectorType& w,const VectorType& d) const
		{
			w.resize(d.size()-1);
			for (size_t i=0;i<w.size();i++)
				w[i] = fabs(d[i+1] - d[i]); 
		}
		
		std::vector<AkimaStruct> akimaStruct_;
		IntervalType interval_;
};

#endif //AKIMA_H_
