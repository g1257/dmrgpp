#ifndef CINC_BASIS_H
#define CINC_BASIS_H
#include "BasisOneSpin.h"
#include "Vector.h"

namespace Dmft {

class Basis {
public:

	Basis(SizeType sites, SizeType nup,SizeType ndown)
	    : nup_(nup),
	      ndown_(ndown),
	      basis1_(sites, nup),
	      basis2_(sites, ndown)
	{}

	SizeType size() const { return basis1_.size()*basis2_.size(); }

private:

	SizeType nup_;
	SizeType ndown_;
	BasisOneSpin basis1_;
	BasisOneSpin basis2_;
};

}
#endif // CINC_BASIS_H
