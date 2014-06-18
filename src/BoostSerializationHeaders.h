#ifndef PSI_BOOST_SER_H_H
#define PSI_BOOST_SER_H_H

#ifdef USE_BOOST
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#else
namespace boost {

namespace serialization {

template<typename T>
class base_object {
};

} // namespace serialization
class archive {

        class text_iarchive {};
};

} // namespace boost
#endif

#endif

