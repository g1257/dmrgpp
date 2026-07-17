#ifndef PSIMAG_KOKKOS_TYPE_H
#define PSIMAG_KOKKOS_TYPE_H

#include <Kokkos_Complex.hpp>

#include <type_traits>

namespace PsimagLite {

// The scalar types that are floating point types and their corresponding std::complex types.
// We need to map std::complex to Kokkos::complex while keeping all the other types which is what
// KokkosType does.
template <typename T> struct KokkosType {
	using type = T;
};

template <typename T> struct KokkosType<std::complex<T>> {
	using type = Kokkos::complex<T>;
};

}

#endif // PSIMAG_KOKKOS_TYPE_H
