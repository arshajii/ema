/// 786

/******************************************************************************/

#pragma once

/******************************************************************************/

#include <map>
#include <array>
#include <vector>
#include <string>
#include <functional>
#include <chrono>
#include <string>
#include <cstdlib>
#include <unordered_map>

#include <sys/types.h>
#include <sys/stat.h>

#include "format.h"

/******************************************************************************/

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

#define prn(f, ...)    fmt::print(f "\n", ##__VA_ARGS__)
#define prnn(...)      fmt::print(__VA_ARGS__)

#define eprn(f, ...)   fmt::print(stderr, f "\n",  ##__VA_ARGS__)
#define eprnn(...)     fmt::print(stderr, __VA_ARGS__)

#define die_if(Q, f, ...)  if (Q) { fmt::print(stderr, f "---exiting!\n",  ##__VA_ARGS__) ; exit(1); }

#ifdef NDEBUG
#define dprn(f, ...)   ;
#define dprnn(...)     ;
#else
#define dprn(f, ...)   { if(getenv("EMADBG")) fmt::print(stderr, f "\n",  ##__VA_ARGS__); }
#define dprnn(...)     { if(getenv("EMADBG")) fmt::print(stderr, __VA_ARGS__); }
#endif

#define cur_time()     chrono::high_resolution_clock::now()
#define elapsed(t)     (chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - (t)).count() / 1000.00)

#define DO(x) for(int _ = 0; _ < (x); _++)
#define DOV(i,x) for(int (i) = 0; (i) < (x); (i)++)

/******************************************************************************/

const int MATE1_TRIM = 7;
const int BC_LEN = 16;

const int ILLUMINA_QUAL_OFFSET = 33;
const int QUAL_BASE = ILLUMINA_QUAL_OFFSET + 1;

const int MIN_READ_SIZE = 32;

const size_t KB = 1024;
const size_t MB = 1024 * KB;
const size_t GB = 1024 * MB;

#define CharToInt(c) ((c) - '0')
#define TwoCharToInt(s) ((10 * CharToInt(s[0])) + CharToInt(s[1])) 
#define PackHaplotagString(s) PackHaplotag(TwoCharToInt(s.substr(1,2)),TwoCharToInt(s.substr(7,2)),TwoCharToInt(s.substr(4,2)),TwoCharToInt(s.substr(10,2)))
#define PackHaplotag(a,b,c,d) ((((uint32_t)(a)) << 24) | (((uint32_t)(c)) << 16) | (((uint32_t)(b)) << 8) | (uint32_t)(d))
#define GenerateAllHaplotagBC(map) for(uint32_t a = 1; a <= 96; ++a) for(uint32_t b = 1; b <= 96; ++b) for(uint32_t c = 1; c <= 96; ++c) for(uint32_t d = 1; d <= 96; ++d) map[PackHaplotag(a,b,c,d)] = 0
#define GenerateAllHaplotagBC_Field(map,f) for(uint32_t a = 1; a <= 96; ++a) for(uint32_t b = 1; b <= 96; ++b) for(uint32_t c = 1; c <= 96; ++c) for(uint32_t d = 1; d <= 96; ++d) map[PackHaplotag(a,b,c,d)].f = 0
/******************************************************************************/

inline char hash_dna(char c)
{
	const static char table[] = { 
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,	0, // 16
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 32
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 64
		0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 
		0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 96
		0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 
		0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  // 128
	};
	return table[c];
}

inline char hash_dna_n(char c)
{
	const static char table[] = { 
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,	0, // 16
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 32
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 64
		0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0, 
		0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 96
		0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 4, 0, 
		0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  // 128
	};
	return table[c];
}


/******************************************************************************/

template<typename K, typename V>
inline size_t estimate_size(const std::map<K, V> &m)
{
	const int overhead = 32;
	return (sizeof(K) + sizeof(V) + overhead) * m.size();
}

template<typename K, typename V>
inline size_t estimate_size(const std::unordered_map<K, V> &M) 
{
	size_t n = M.bucket_count();
	float m = M.max_load_factor();
	if (m > 1.0) n *= m;
	return (M.size() * sizeof(V) + n * (sizeof(size_t) + sizeof(void*)));
}

inline int stat_dir(const std::string &path)
{
	struct stat path_stat;
	int s = stat(path.c_str(), &path_stat);
	if (s != 0) {
		return 0;
	}
	if (S_ISDIR(path_stat.st_mode)) {
		return 1;
	}
	if (S_ISREG(path_stat.st_mode)) {
		return 2;
	} 
	return 3;
}
