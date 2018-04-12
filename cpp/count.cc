/// 786

#include <unordered_map>
#include <string>
#include <fstream>
#include <iostream>

#include "common.h"

using namespace std;



/******************************************************************************/

void dump_map(map<string, int64_t> &full_counts, FILE *fo)
{
	static int count = 0;

	// eprn("current map size: {:n}", estimate_size(full_counts));

	int64_t total_buckets = full_counts.size();
	size_t c = fwrite(&total_buckets, 8, 1, fo);
	die_if(c != 1, "fwrite failed");
	for (auto &bcd: full_counts) {
		c = fwrite(bcd.first.data(), 1, BC_LEN, fo);
		die_if(c != BC_LEN, "fwrite failed");
		c = fwrite(&bcd.second, 8, 1, fo);
		die_if(c != 1, "fwrite failed");
	}
	fflush(fo);
	full_counts.clear();
	eprn(":: Dumped block {:n}", ++count);
}

/******************************************************************************/

EXTERNC void count(
	const char *known_barcodes_path,
	const char *output_prefix,
	const size_t max_map_size)
{
	auto TOT = cur_time();

	unordered_map<uint32_t, int64_t> counts;
	map<string, int64_t> full_counts;

	string s, q;

	auto T = cur_time();
	ifstream fin(known_barcodes_path);
	die_if(fin.fail(), "Cannot open file {}", known_barcodes_path);
	uint32_t barcode;
	while (getline(fin, s)) {
		barcode = 0;
		DO(BC_LEN) barcode = (barcode << 2) | hash_dna(s[_]);
		die_if(barcode == 0, "Invalid barcode AAA...AA whitelisted");
		counts[barcode] = 0;
	}
	fin.close();
	eprn(":: Loading 10X took {:.1f} s", elapsed(T)); T = cur_time();

	FILE *f_full = fopen(fmt::format("{}.ema-fcnt", output_prefix).c_str(), "wb");
	die_if(f_full == NULL, "Cannot open file {}", fmt::format("{}.ema-fcnt", output_prefix));
	FILE *f_nice = fopen(fmt::format("{}.ema-ncnt", output_prefix).c_str(), "wb");
	die_if(f_nice == NULL, "Cannot open file {}", fmt::format("{}.ema-ncnt", output_prefix));

	string b = string(BC_LEN, '#');
	int64_t total_reads = 0;
	int64_t nice_reads = 0;
	int64_t ignored_reads = 0;
	int64_t sz = 0;
	while (getline(cin, s)) { // N R1 Q1 R2 Q2
		sz += s.size() + 1;
		getline(cin, s); sz += s.size() + 1;
		getline(cin, q); sz += q.size() + 1;
		getline(cin, q); sz += q.size() + 1;

		bool process = s.size() >= MIN_READ_SIZE;

		bool has_n = 0;
		barcode = 0;
		if (process) DO(BC_LEN) { // 0..33 34..
			//die_if((s[_] == 'N' && q[_] != '#') || (s[_] != 'N' && q[_] == '#'), 
			//	"# quality score does not match N ({} vs. {})", s[_], q[_]);
			if (q[_] < ILLUMINA_QUAL_OFFSET) {
				process = false;
				eprn("Ignoring long read--- quality score {} less than {}", q, ILLUMINA_QUAL_OFFSET);
				break;
			}
			if (q[_] - ILLUMINA_QUAL_OFFSET >= QUAL_BASE) {
				//eprn("Trimming quality score {} to {}", q[_], char(ILLUMINA_QUAL_OFFSET + QUAL_BASE - 1));
				q[_] = ILLUMINA_QUAL_OFFSET + QUAL_BASE - 1;
			}

			b[_] = hash_dna_n(s[_]) * QUAL_BASE + min(QUAL_BASE - 1, q[_] - ILLUMINA_QUAL_OFFSET);
			barcode = (barcode << 2) | hash_dna(s[_]);
			has_n |= (s[_] == 'N');
		}
		if (process) {
			if (!has_n) {
				auto it = counts.find(barcode);
				if (it != counts.end()) {
					it->second++;
					nice_reads++;
				}
			}
			int cnt = full_counts[b]++;
			if (!cnt && estimate_size(full_counts) >= max_map_size) { // new element
				dump_map(full_counts, f_full);
			}
			total_reads++;
		} else {
			ignored_reads++;
		}

		DO(4) { getline(cin, s); sz += s.size() + 1; }
	}
	eprn(":: Counting took {:.1f} s", elapsed(T)); T = cur_time();
	eprn(":: Reads with OK barcode: {:n} out of {:n}", nice_reads, total_reads);
	eprn(":: Ignored {:n} reads", ignored_reads);

	int64_t tenx_buckets = 0;
	for (auto &bcd: counts) if (bcd.second) {
		tenx_buckets++;
	}
	size_t c = fwrite(&tenx_buckets, 8, 1, f_nice);
	die_if(c != 1, "fwrite failed");
	for (auto &bcd: counts) if (bcd.second) {
		c = fwrite(&bcd.first, 4, 1, f_nice);
		die_if(c != 1, "fwrite failed");
		c = fwrite(&bcd.second, 8, 1, f_nice);
		die_if(c != 1, "fwrite failed");
	}
	fclose(f_nice);

	dump_map(full_counts, f_full);
	fclose(f_full);
	eprn(":: Printing took {:.1f} s", elapsed(T)); T = cur_time();

	eprn(":: Processed {:n} reads ({:n} MB uncompressed) in {:n} s",
		total_reads,
		sz / (1024 * 1024), (int)elapsed(TOT));
}
