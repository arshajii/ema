/// 786

#include <unordered_map>
#include <map>
#include <string>
#include <thread>
#include <queue>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <valarray>
#include <mutex>

#include <fcntl.h>
#include <unistd.h>

#include "common.h"

using namespace std;

/******************************************************************************/

const double BC_CONF_THRESH = 0.975;

const int NOCHANGE = 0;
const int H1CHANGE = 1;
const int H2CHANGE = 2;
const int NOBUCKET = 3;

/******************************************************************************/

struct Count {
	int64_t n_reads; // How many reads falling in here?
	double prior;
	int bucket; // Which bucket
	
	Count(): n_reads(0), bucket(0), prior(0) {}
};

struct Buffer {
	char *buff;
	size_t size;
};

mutex correct_mutex;

/******************************************************************************/

static double _phred_probs[128];

void initialize_probs()
{
	for (int i = 0; i < 128; i++) {
		::_phred_probs[i] = pow(10.0, - min(QUAL_BASE - 1, i) / 10.0);
	}
}

double get_prob(char c)
{
	return ::_phred_probs[c];
}

/******************************************************************************/

void correct_barcode(int thread, int max_threads,
	unordered_map<uint32_t, Count> &known_counts,
	map<string, pair<uint32_t, int64_t>> &full_counts,
	map<string, uint32_t> &final_corrected_counts,
	vector<int64_t> &final_stats,
	bool do_h2)
{
	auto T = cur_time();

	map<string, uint32_t> corrected;
	vector<int64_t> stats(4, 0);

	int64_t it_c = 0;
	for (auto &it: full_counts) {
		/// TODO: better threading with std::map if it's possible at all
		if (it_c++ % max_threads != thread) { 
			continue;
		}

		const string &q = it.first;
		int64_t fc = it.second.second;

		// Is it already corrected?
		// auto itt = corrected_counts.find(q);
		// if (itt != corrected_counts.end()) {
		// 	it.second.first = itt->second.first;
		// 	stats[itt->second.second] += fc;
		// 	continue;
		// }

		
		uint32_t barcode = 0;
		it.second.first = 0; // Assign empty barcode
		int ns = 0; // how many Ns?
		DO(BC_LEN) {
			barcode = (barcode << 2) | (uint8_t(q[_]) / QUAL_BASE);
			ns += (uint8_t(q[_]) % QUAL_BASE == 0);
		}
		if (ns > 1) {
			stats[NOBUCKET] += fc;
			continue;
		}

		auto cit = (ns == 0) ? known_counts.find(barcode) : known_counts.end();
		uint32_t max_barcode = 0;
		double max_p = -1;
		double total = 0;

		int type = NOBUCKET;
		if (cit != known_counts.end()) {
			max_p = cit->second.prior;
			max_barcode = barcode;
			total += max_p;
			type = NOCHANGE;
			if (do_h2) DOV(i1, BC_LEN) DOV(j1, 4) {
				if (j1 == uint8_t(q[i1]) / QUAL_BASE)
					continue;
				for (int i2 = i1 + 1; i2 < BC_LEN; i2++) DOV(j2, 4) {
					if (j2 == uint8_t(q[i2]) / QUAL_BASE)
						continue;
					uint32_t bcd = barcode 
						& ~(3 << ((BC_LEN - i1 - 1) * 2)) 
						& ~(3 << ((BC_LEN - i2 - 1) * 2)) 
						& ((1ll << (BC_LEN * 2)) - 1);
					bcd |= (j1 << ((BC_LEN - i1 - 1) * 2));
					bcd |= (j2 << ((BC_LEN - i2 - 1) * 2));
					cit = known_counts.find(bcd);
					if (cit != known_counts.end()) {
						double p1 = get_prob(max(3.0, uint8_t(q[i1]) % QUAL_BASE - 1.0)),
						       p2 = get_prob(max(3.0, uint8_t(q[i2]) % QUAL_BASE - 1.0));
						double p = cit->second.prior * (p1 * p2);
						total += p;
						if (p > max_p) {
							max_p = p;
							max_barcode = bcd;
							type = H2CHANGE;
						}							
					}
				}
			}
		} else if (ns <= 1) {
			DOV(i, BC_LEN) {
				if (ns && uint8_t(q[i]) % QUAL_BASE != 0) 
					continue;
				DOV(j, 4) {
					if (ns == 0 && j == uint8_t(q[i]) / QUAL_BASE)
						continue;
					uint32_t bcd = barcode & ~(3 << ((BC_LEN - i - 1) * 2)) & ((1ll << (BC_LEN * 2)) - 1);
					bcd |= (j << ((BC_LEN - i - 1) * 2));
					
					auto cit = known_counts.find(bcd);
					if (cit != known_counts.end()) {
						double p = cit->second.prior * get_prob(uint8_t(q[i]) % QUAL_BASE);
						total += p;
						if (p > max_p) {
							max_p = p;
							max_barcode = bcd;
							type = H1CHANGE;
						}							
					}
				}
			}
		} 

		if (max_p / total > BC_CONF_THRESH) {
			it.second.first = max_barcode;
			if (type == H1CHANGE || type == H2CHANGE) {
				corrected[q] = max_barcode;
			}
		} else {
			type = NOBUCKET;
		}

		stats[type] += fc;
	}

	{
		lock_guard<mutex> guard(correct_mutex);

		for (int i = 0; i < 4; i++) {
			final_stats[i] += stats[i];
		}
		for (auto &x: corrected) {
			final_corrected_counts[x.first] = x.second;
		}
		it_c = 0;
		for (auto &c: full_counts) {
			if (it_c++ % max_threads != thread) { 
				continue;
			}
			auto new_bcd = c.second.first;
			if (new_bcd != 0) {
				// assert(known_counts.find(new_bcd) != known_counts.end());
				known_counts[new_bcd].n_reads += c.second.second;
			}
		}
	}

	eprn("  :: Thread {} took {} s", thread, elapsed(T));
}

/******************************************************************************/

void load_barcode_count(const string &path, unordered_map<uint32_t, Count> &counts)
{
	eprn(":: Loading known counts file {} ... ", path);
	FILE *f = fopen(path.c_str(), "rb");
	assert(f != NULL);
	int64_t total;
	size_t c = fread(&total, 8, 1, f);
	assert(c == 1);
	while (total--) {
		uint32_t bcd;
		c = fread(&bcd, 4, 1, f);
		assert(c == 1);
		int64_t cnt;
		c = fread(&cnt, 8, 1, f);
		assert(c == 1);
		counts[bcd].prior += cnt;
	}
	fclose(f);
}

void load_and_correct_full_count(const string &path, 
	unordered_map<uint32_t, Count> &known_counts,
	map<string, uint32_t> &corrected_counts,
	vector<int64_t> &stats,
	bool do_h2,
	const int nthreads)
{
	eprn(":: Loading full counts file {} ... ", path);

	FILE *f = fopen(path.c_str(), "rb");
	assert(f != NULL);

	int64_t total, dump = 0;
	size_t c;
	string b = string(BC_LEN, '#');
	while (fread(&total, 8, 1, f) == 1) {
		auto T = cur_time();
		
		map<string, pair<uint32_t, int64_t>> full_counts;
		while (total--) {
			c = fread(&b[0], 1, BC_LEN, f);
			assert(c == BC_LEN);
			int64_t cnt;
			c = fread(&cnt, 8, 1, f);
			assert(c == 1);
			full_counts[b].first = 0;
			full_counts[b].second += cnt;
		}
		eprn("::: Loading dump {} done in {} s", ++dump, elapsed(T)); T = cur_time();

		vector<thread> threads;

		for (int i = 0; i < nthreads; i++) {
			threads.push_back(thread(
				correct_barcode, 
				i, nthreads,
				ref(known_counts), 
				ref(full_counts), 
				ref(corrected_counts),
				ref(stats),
				do_h2
			));
		}
		for (auto &t: threads) {
			t.join(); 
		}
		eprn("::: Correcting done in {} s", elapsed(T)); T = cur_time();
		eprn("--> corrected counts size: {:n}", estimate_size(corrected_counts));
	}
	fclose(f);
}


/******************************************************************************/

EXTERNC void correct(
	const char *known_barcodes_path, 
	const char **input_prefix, 
	const int input_prefix_size,
	const char *output_dir, 
	const char do_h2,
	const size_t buffer_size,
	const int nthreads, 
	const int nbuckets) 
{
	auto T = cur_time();
	initialize_probs();
	eprn(":: Bucketing {} inputs into {} files with {} threads", 
		input_prefix_size, nbuckets, nthreads);
	
/// 1. Load known counts
	unordered_map<uint32_t, Count> known_counts;
	string s, q, n, r;
	ifstream fin(known_barcodes_path);
	assert(!fin.fail());
	while (getline(fin, s)) {
		uint32_t barcode = 0;
		DO(BC_LEN) barcode = (barcode << 2) | hash_dna(s[_]);
		assert(barcode != 0);
		known_counts[barcode].prior = 0;
	}
	fin.close();

	for (int fi = 0; fi < input_prefix_size; fi++) {
		load_barcode_count(fmt::format("{}.ema-ncnt", input_prefix[fi]), known_counts);
	}
	double total_counts = 0;
	for (auto &c: known_counts) {
		total_counts += c.second.prior + 1;
	}
	for (auto &c: known_counts) {
		c.second.prior = (c.second.prior + 1) / total_counts;
	}
	eprn(":: Loading known counts ... done in {:.1f} s", elapsed(T)); T = cur_time();

/// 2. Load and correct full counts
	map<string, uint32_t> corrected_counts;
	vector<int64_t> stats(4, 0);
	for (int fi = 0; fi < input_prefix_size; fi++) {
		load_and_correct_full_count(
			fmt::format("{}.ema-fcnt", input_prefix[fi]), 
			known_counts,
			corrected_counts, 
			stats, 
			do_h2,
			nthreads
		);
	}
	eprn(":: Stats: no change: {:n} \n"
		 "         no barcode: {:n} \n"
		 "       H1-corrected: {:n} \n"
		 "       H2-corrected: {:n} ", 
		 stats[NOCHANGE], stats[NOBUCKET], stats[H1CHANGE], stats[H2CHANGE]);
	eprn(":: Corrected map size: {:n} MB", estimate_size(corrected_counts) / MB);
	eprn(":: Corrected size: {:n}", corrected_counts.size());
	eprn(":: Correcting barcodes ... done in {:.1f} s", elapsed(T)); T = cur_time();

/// 3. Set up bucket boundaries
	// 3a. Create file handles and buffers
	assert(S_ISDIR(stat_file(output_dir)));
	struct PQFile {
		int64_t size; // current size
		FILE *file;
		char *buf;
		size_t buf_size;
	};
	vector<PQFile> files { PQFile { // Non-bucketed reads
		.size = 0,
		.file = fopen(fmt::format("{}/ema-bin-nobc", output_dir).c_str(), "wb"),
		.buf = (char*)malloc(buffer_size + 10 * KB),
		.buf_size = 0
	}};
	assert(files.back().buf != NULL);
	assert(files.back().file != NULL);

	// Have the smallest file at the top
	auto PQComp = [&](int a, int b) { 
		return tie(files[a].size, a) > tie(files[b].size, b); 
	}; 
	priority_queue<int, vector<int>, decltype(PQComp)> pq(PQComp);
	for (int i = 0; i < nbuckets; i++) {
		files.push_back({
			.size = 0,
			.file = fopen(fmt::format("{}/ema-bin-{:03}", output_dir, i).c_str(), "wb"),
			.buf = (char*)malloc(buffer_size + 10 * KB),
			.buf_size = 0
		});
		assert(files.back().buf != NULL);
		assert(files.back().file != NULL);
		pq.push(files.size() - 1);
	}
	/// 3b. Bucket the barcodes
	for (auto &cnt: known_counts) {
		int fidx = pq.top(); pq.pop();
		files[fidx].size += cnt.second.n_reads;
		cnt.second.bucket = fidx;
		pq.push(fidx);
	}
	eprn(":: File assignment ... done in {:.1f} s", elapsed(T)); T = cur_time();

	for (int i = 0; i < files.size(); i++) {
		eprn("{} -> {:n}", i, files[i].size);
	}
	// exit(0);
	
/// 4. Write buckets
	char bcd[BC_LEN + 1];
	bcd[BC_LEN] = 0;
	string b = string(BC_LEN, '#');

	// ifstream cin("../../ema/data/inter.fq");
	while (getline(cin, n)) { 
		getline(cin, r);
		getline(cin, q);
		getline(cin, q);

		uint32_t barcode;
		bool has_n = 0;
		DO(BC_LEN) {
			barcode = (barcode << 2) | hash_dna(r[_]);
			has_n |= (r[_] == 'N');
			b[_] = hash_dna(r[_]) * QUAL_BASE + min(QUAL_BASE - 1, q[_] - ILLUMINA_QUAL_OFFSET);
		}
		
		int fidx;
		auto cit = corrected_counts.find(b);
		if (cit != corrected_counts.end()) {
			barcode = cit->second;
			has_n = 0;
		} 
		auto kit = has_n ? known_counts.end() : known_counts.find(barcode);
		if (kit != known_counts.end()) {
			fidx = kit->second.bucket;
		} else {
			barcode = 0;
			fidx = 0;
		}
		
		char *buff = files[fidx].buf;
		size_t &buffi = files[fidx].buf_size;

		// Barcode
		if (barcode != 0) {
			DO(BC_LEN) bcd[BC_LEN - _ - 1] = "ACGT"[barcode & 3], barcode >>= 2;
			memcpy(buff + buffi, bcd, BC_LEN); 
			buffi += BC_LEN;
			buff[buffi++] = ' ';
		}

		// Name
		for (auto c: n) {
			if (isspace(c)) break; 
			buff[buffi++] = c;
		}
		buff[buffi++] = ' ';

		// Trimmed read
		memcpy(buff + buffi, r.c_str() + BC_LEN + MATE1_TRIM, r.size() - BC_LEN - MATE1_TRIM); 
		buffi += r.size() - BC_LEN - MATE1_TRIM;	
		buff[buffi++] = ' ';

		// Trimmed quality
		memcpy(buff + buffi, q.c_str() + BC_LEN + MATE1_TRIM, q.size() - BC_LEN - MATE1_TRIM); 
		buffi += r.size() - BC_LEN - MATE1_TRIM;
		buff[buffi++] = ' ';
		
		getline(cin, s); 
		getline(cin, s); 
		// Pair read
		memcpy(buff + buffi, s.c_str(), s.size()); 
		buffi += s.size();
		buff[buffi++] = ' ';
		getline(cin, s); 
		getline(cin, s); 
		// Pair quality
		memcpy(buff + buffi, s.c_str(), s.size()); 
		buffi += s.size();
		buff[buffi++] = '\n';

		if (buffi >= buffer_size) {
			fwrite(buff, 1, buffi, files[fidx].file);
			buffi = 0;
		}
	}
	for (int fidx = 0; fidx < files.size(); fidx++) {
		fwrite(files[fidx].buf, 1, files[fidx].buf_size, files[fidx].file);
		fclose(files[fidx].file);
		free(files[fidx].buf);
	}
	eprn("Writing barcodes ... done in {:.1f} s", elapsed(T)); T = cur_time();
}

