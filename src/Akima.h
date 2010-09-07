#ifndef AKIMA_H_
#define AKIMA_H_

//! A class to interpolate using akima spline
class Akima {
	struct AkimaStruct {
		RealType x0,x1,a0,a1,a2,a3;
	};

	public:
		Akima(const VectorType& x,const VectorType& s)
		{
			std::vector<RealType> sprime;
			calculateSprime(sprime,x,s);
			for (size_t i=0;i<n;i++) {
				RealType u = x[i+1]-x[i];
				RealType a2 = u*(sprime[i+1]-sprime[i])-2*(s[i+1]-s[i]));
				RealType a3 = 3*u*(s[i+1]-s[i])-u*u*(sprime[i+1]+2*sprime[i]);
				RealType ucube = u*u*u;
				AkimaStruct ak(x[i],x[i+1],s[i],sprime[i],a2/ucube,a3/ucube);
				akimaStruct_.push_back(ak);
			}
		}
			
		RealType operator()(const RealType& x) const
		{
			size_t i =findRange(x);
			return functionFor(x,i);
		}

	private:
		
		RealType functionFor(const RealType& x,size_t i) const
		{
			const AkimaStruct& ak = akimaStruct_[i];
			RealType u = x - ak.x0;
			return ak.a0+ak.a1*u+ak.a2*u*u+ak.a3*u*u*u;
		}
		
		size_t findRange(const RealType& x) const
		{
			if (x<akimaStruct_[0].x0 || x>akimaStruct_[0].x1)
				throw std::runtime_error("Akima: out of range\n");
			size_t l = akimaStruct_.size();
			size_t i = l/2;
			while (x<akimaStruct_[i].x0) {
				i = i/2;
			}
			while(x>akimaStruct_[i].x1) {
				i = (l+i)/2;
			}
			return i;
		}
		
		std::vector<AkimaStruct> akimaStruct_;
};

#endif //AKIMA_H_
