/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
/** \ingroup PsimagLite */
/*@{*/

/*!
 *
 *
 */
#ifndef MAP_HEADER_H
#define MAP_HEADER_H
#include <map>
#include "AllocatorCpu.h"

namespace PsimagLite {

template<typename Key,typename T,typename Compare = std::less<Key> >
class Map {
public:

	typedef std::map<Key,T,Compare,typename Allocator<std::pair<const Key,T> >::Type> Type;
}; // class Map

template<typename T>
class IsMapLike {
public:
	enum {True = false};
};

template<typename Key,typename T,typename Compare>
class IsMapLike<std::map<Key,T,Compare,typename Allocator<std::pair<const Key,T> >::Type> > {
public:
	enum {True = true};
};

template<typename MapType>
typename HasType<IsMapLike<MapType>::True,void>::Type
printMap(std::ostream& os, const MapType& x,const String& label)
{
	typedef typename MapType::const_iterator MapIteratorType;
	for (MapIteratorType it = x.begin();it!=x.end();++it) {
		os<<label<<"["<<it->first<<"]="<<it->second<<"\n";
	}
}

} // namespace PsimagLite

/*@}*/
#endif // MAP_HEADER_H
