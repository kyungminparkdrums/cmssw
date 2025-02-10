#ifndef DATAFORMATS_L1TPARTICLEFLOW_ENCODING_H
#define DATAFORMATS_L1TPARTICLEFLOW_ENCODING_H

#include <cassert>
#include "DataFormats/L1TParticleFlow/interface/datatypes.h"

template <typename U, typename T>
inline void pack_into_bits(U& u, unsigned int& start, const T& data) {
  const unsigned int w = T::width;
  u(start + w - 1, start) = data(w - 1, 0);
  start += w;
}

template <typename U, typename T>
inline void unpack_from_bits(const U& u, unsigned int& start, T& data) {
  const unsigned int w = T::width;
  data(w - 1, 0) = u(start + w - 1, start);
  start += w;
}

template <typename U>
inline void pack_bool_into_bits(U& u, unsigned int& start, bool data) {
  u[start++] = data;
}

template <typename U>
inline void unpack_bool_from_bits(const U& u, unsigned int& start, bool& data) {
  data = u[start++];
}

enum class PackingStrategy { DEFAULT, BARREL, ENDCAP };

// Helper function that calls the correct unpack method based on user input
template <typename T, int NB, PackingStrategy METHOD = PackingStrategy::DEFAULT>
inline auto unpack_helper(const ap_uint<NB>& data) {
  if constexpr (METHOD == PackingStrategy::BARREL) {
    static_assert(T::BITWIDTH_BARREL <= NB, "NB Type is too small for the object");
    return T::unpack_barrel(data);
  } else if constexpr (METHOD == PackingStrategy::ENDCAP) {
    static_assert(T::BITWIDTH_ENDCAP <= NB, "NB Type is too small for the object");
    return T::unpack_endcap(data);
  } else {
    static_assert(T::BITWIDTH <= NB, "NB Type is too small for the object");
    return T::unpack(data);
  }
}

template <typename T, int NB, PackingStrategy METHOD = PackingStrategy::DEFAULT>
inline auto pack_helper(const T& obj) {
  if constexpr (METHOD == PackingStrategy::BARREL) {
    static_assert(T::BITWIDTH_BARREL <= NB, "NB Type is too small for the object");

    return obj.pack_barrel();
  } else if constexpr (METHOD == PackingStrategy::ENDCAP) {
    static_assert(T::BITWIDTH_ENDCAP <= NB, "NB Type is too small for the object");

    return obj.pack_endcap();
  } else {
    static_assert(T::BITWIDTH <= NB, "NB Type is too small for the object");

    return obj.pack();
  }
}

template <typename T, int NB, PackingStrategy METHOD>
inline auto pack_slim_helper(const T& obj) {
  if constexpr (METHOD == PackingStrategy::BARREL) {
    static_assert(T::BITWIDTH_BARREL_SLIM <= NB, "NB Type is too small for the object");
    return obj.pack_barrel_slim();
  } else if constexpr (METHOD == PackingStrategy::ENDCAP) {
    static_assert(T::BITWIDTH_ENDCAP_SLIM <= NB, "NB Type is too small for the object");
    return obj.pack_endcap_slim();
  } else {
    static_assert(T::BITWIDTH_SLIM <= NB, "NB Type is too small for the object");
    return obj.pack_slim();
  }
}

template <unsigned int N, PackingStrategy METHOD = PackingStrategy::DEFAULT, unsigned int OFFS = 0, typename T, int NB>
inline void l1pf_pattern_pack(const T objs[N], ap_uint<NB> data[]) {
#ifdef __SYNTHESIS__
#pragma HLS inline
#pragma HLS inline region recursive
#endif
  for (unsigned int i = 0; i < N; ++i) {
#ifdef __SYNTHESIS__
#pragma HLS unroll
#endif
    data[i + OFFS] = pack_helper<T, NB, METHOD>(objs[i]);
  }
}

template <unsigned int N, PackingStrategy METHOD = PackingStrategy::DEFAULT, unsigned int OFFS = 0, typename T, int NB>
inline void l1pf_pattern_unpack(const ap_uint<NB> data[], T objs[N]) {
#ifdef __SYNTHESIS__
#pragma HLS inline
#pragma HLS inline region recursive
#endif
  for (unsigned int i = 0; i < N; ++i) {
#ifdef __SYNTHESIS__
#pragma HLS unroll
#endif
    objs[i] = unpack_helper<T, NB, METHOD>(data[i + OFFS]);
  }
}

template <unsigned int N, PackingStrategy METHOD = PackingStrategy::DEFAULT, unsigned int OFFS = 0, typename T, int NB>
inline void l1pf_pattern_pack_slim(const T objs[N], ap_uint<NB> data[]) {
#ifdef __SYNTHESIS__
#pragma HLS inline
#pragma HLS inline region recursive
#endif
  for (unsigned int i = 0; i < N; ++i) {
#ifdef __SYNTHESIS__
#pragma HLS unroll
#endif
    data[i + OFFS] = pack_slim_helper<T, NB, METHOD>(objs[i]);
  }
}

template <unsigned int N, PackingStrategy METHOD = PackingStrategy::DEFAULT, unsigned int OFFS = 0, typename T, int NB>
inline void l1pf_pattern_unpack_slim(const ap_uint<NB> data[], T objs[N]) {
#ifdef __SYNTHESIS__
#pragma HLS inline
#pragma HLS inline region recursive
#endif
  for (unsigned int i = 0; i < N; ++i) {
#ifdef __SYNTHESIS__
#pragma HLS unroll
#endif
    objs[i] = unpack_helper<T, NB, METHOD>(data[i + OFFS]);
  }
}

#endif
