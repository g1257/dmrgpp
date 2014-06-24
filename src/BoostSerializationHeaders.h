#ifndef PSI_BOOST_SER_H_H
#define PSI_BOOST_SER_H_H
#include <exception>

#ifdef USE_BOOST
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_free.hpp>
#else
#include <string>
#include <stdexcept>
namespace boost {

void noboost(std::string msg)
{
	std::string str("ERROR: Please compile with -DUSE_BOOST for ");
	str += msg + " to work\n";
	throw std::runtime_error(str);
}

namespace serialization {

template<typename T>
class base_object {
};

} // namespace serialization

namespace archive {

class text_iarchive {

public:

	text_iarchive(std::ifstream&)
	{
		boost::noboost("text_oarchive");
	}
};

template<typename T>
text_iarchive& operator>>(text_iarchive& ia, const T& t)
{
	boost::noboost("operator>> text_iarchive");
	return ia;
}

class text_oarchive {

public:

	text_oarchive(std::ofstream&)
	{
		boost::noboost("text_oarchive");
	}

};

template<typename T>
text_oarchive& operator<<(text_oarchive& oa, const T& t)
{
	boost::noboost("operator<< text_oarchive");
	return oa;
}

} // namespace archive

} // namespace boost
#endif

#endif

