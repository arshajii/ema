/// 786

#include <unordered_map>
#include <string>
#include <fstream>
#include <iostream>

#include "common.h"

using namespace std;

/******************************************************************************/

void count(const string &path) 
{	
	unordered_map<uint32_t, int64_t> counts;
	map<string, int64_t> full_counts;
	
	string s, q;

	auto T = cur_time();
	auto TOT = cur_time();
	ifstream fin(path.c_str());
	uint32_t barcode;
	while (getline(fin, s)) {
		barcode = 0;
		DO(BC_LEN) barcode = (barcode << 2) | hash_dna(s[_]);
		assert(barcode != 0);
		counts[barcode] = 0;
	}
	fin.close();
	eprn(":: Loading 10X took {:.1f} s", elapsed(T)); T = cur_time();

	string b = string(BC_LEN, '#');
	int64_t total_reads = 0;
	int64_t nice_reads = 0;
	int64_t sz = 0;
	while (getline(cin, s)) { // N R1 Q1 R2 Q2		
		sz += s.size() + 1;
		getline(cin, s); sz += s.size() + 1; 
		getline(cin, q); sz += q.size() + 1; 
		getline(cin, q); sz += q.size() + 1; 
		bool has_n = 0;
		barcode = 0;
		DO(BC_LEN) { // 0..33 34..
			if (s[_] == 'N') assert(q[_] == '#');
			if (q[_] == '#') assert(s[_] == 'N');
			assert(q[_] >= ILLUMINA_QUAL_OFFSET);
			b[_] = hash_dna(s[_]) * QUAL_BASE + min(QUAL_BASE - 1, q[_] - ILLUMINA_QUAL_OFFSET);
			barcode = (barcode << 2) | hash_dna(s[_]);
			has_n |= (s[_] == 'N');
		}
		// bool not_added = 1;
		if (!has_n) {
			auto it = counts.find(barcode);
			if (it != counts.end()) {
				it->second++;
				nice_reads++;
				// not_added = 0;
			}
		}
		// if (not_added) {
		full_counts[b]++;
		// }
		getline(cin, s); sz += s.size() + 1;
		getline(cin, s); sz += s.size() + 1; 
		getline(cin, s); sz += s.size() + 1;
		getline(cin, s); sz += s.size() + 1; 
		total_reads++;
	}
	eprn(":: Counting took {:.1f} s", elapsed(T)); T = cur_time();
	eprn(":: 10x-compatible reads: {:n} out of {:n}", nice_reads, total_reads);
	eprn(":: Total buckets: 10x {:n}, total+qualities {:n}", counts.size(), full_counts.size()); 

	FILE *fo = stdout;
	fwrite(&total_reads, 8, 1, fo);
	int64_t tenx_buckets = 0;
	for (auto &bcd: counts) if (bcd.second) 
		tenx_buckets++;
	fwrite(&tenx_buckets, 8, 1, fo);
	for (auto &bcd: counts) if (bcd.second) {
		fwrite(&bcd.first, 4, 1, fo);
		fwrite(&bcd.second, 8, 1, fo);
	}
	int64_t total_buckets = full_counts.size();
	fwrite(&total_buckets, 8, 1, fo);
	for (auto &bcd: full_counts) {
		fwrite(bcd.first.data(), 1, BC_LEN, fo);
		fwrite(&bcd.second, 8, 1, fo);
	}
	fflush(fo);
	eprn(":: Printing took {:.1f} s", elapsed(T)); T = cur_time();
	eprn(":: Processed {:n} reads ({:n} MB) in {:n} s", total_reads, 
		sz / (1024 * 1024), (int)elapsed(TOT));
}
