/// 786

#include <unordered_map>
#include <map>
#include <string>
#include <thread>
#include <queue>
#include <iostream>
#include <algorithm>
#include <fstream>

#include <glob.h>
#include <fcntl.h>
#include <unistd.h>

#include "common.h"

using namespace std;

/******************************************************************************/

struct full_count_t {
	int64_t count;
	uint32_t real_bcd; // 0 is no barcode! requires AAA.AAA is not in the main list
	full_count_t(): count(0), real_bcd(0) {}
};

struct stats_t {
	int64_t total_passed, not_bucketed, corrected_h1, corrected_h2;
	int64_t total_passed_bcd, not_bucketed_bcd, corrected_h1_bcd, corrected_h2_bcd;
};

/******************************************************************************/

void correct_barcode(int thread, int max_threads,
	const unordered_map<uint32_t, double> &counts,
	map<string, full_count_t> &full_counts,
	stats_t &s,
	bool do_h2)
{
	const double BC_CONF_THRESH = 0.975;
	auto T = cur_time();

	s = stats_t {0, 0, 0, 0, 0, 0, 0, 0};

	int64_t it_c = 0;
	for (auto &it: full_counts) {
		it_c++;
		if (it_c % max_threads != thread)
			continue;

		const string &q = it.first;
		full_count_t &fc = it.second;
		s.total_passed += fc.count;
		s.total_passed_bcd++;

		uint32_t barcode = 0;
		int ns = 0;
		DO(BC_LEN) {
			// assert(uint8_t(q[_]) / QUAL_BASE < 4);
			barcode = (barcode << 2) | (uint8_t(q[_]) / QUAL_BASE);
			ns += (uint8_t(q[_]) % QUAL_BASE == 0);
		}
		if (ns > 1) {
			s.not_bucketed += fc.count;
			s.not_bucketed_bcd++;
			fc.real_bcd = 0;
			continue; // Nothing to find here!
		}

		auto cit = (ns == 0) ? counts.find(barcode) : counts.end();
		uint32_t max_barcode = 0;
		double max_p = -1; // check
		double total = 0;

		int cor_type = 3;
		if (cit != counts.end()) {
			max_p = cit->second;
			max_barcode = barcode;
			total += max_p;
			cor_type = 0;
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
					cit = counts.find(bcd);
					if (cit != counts.end()) {
						double p1 = get_prob(max(3.0, uint8_t(q[i1]) % QUAL_BASE - 1.0)),
						       p2 = get_prob(max(3.0, uint8_t(q[i2]) % QUAL_BASE - 1.0));
						double p = cit->second * (p1 * p2);
						total += p;
						if (p > max_p) {
							max_p = p;
							max_barcode = bcd;
							cor_type = 2;
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
					

					auto cit = counts.find(bcd);
					if (cit != counts.end()) {
						double p = cit->second * get_prob(uint8_t(q[i]) % QUAL_BASE);
						total += p;
						if (p > max_p) {
							// uint32_t z;
							// z= barcode; string x; DO(BC_LEN)x="ACGT"[z&3]+x,z>>=2;
							// z= bcd; string y; DO(BC_LEN)y="ACGT"[z&3]+y,z>>=2;
							// // eprn("{} ({}->{}) --> {}", x, i, "ACGT"[j], y);
							// x[i]="ACGT"[j];
							// assert(x==y);

							max_p = p;
							max_barcode = bcd;
							cor_type = 1;
						}							
					}
				}
			}
		} 
		// if (cor_type == 1)
			// eprn("came here, {}, p={}, total={}", cor_type, max_p, total);

		if (max_p / total > BC_CONF_THRESH) {
			fc.real_bcd = max_barcode;
			// if (cor_type == 1) {
			// 	auto z= barcode; string x; DO(BC_LEN)x="ACGT"[z&3]+x,z>>=2;
			// 	z= max_barcode; string y; DO(BC_LEN)y="ACGT"[z&3]+y,z>>=2;
			// 	printf("%s -> %s %lf ~~ %lf\n", x.c_str(),y.c_str(), max_p, max_p/total);
			// }
		} else {
			fc.real_bcd = 0;
			cor_type = 3;
		}

		if (cor_type == 1) s.corrected_h1 += fc.count, s.corrected_h1_bcd++;
		if (cor_type == 2) s.corrected_h2 += fc.count, s.corrected_h2_bcd++; 
		if (cor_type == 3) s.not_bucketed += fc.count, s.not_bucketed_bcd++; 
	}

	eprn("  :: Thread {} took {} s", thread, elapsed(T));
}

/******************************************************************************/

void correct(const string &path, 
	const string &prefix, 
	const int nthreads, 
	const int nfiles) 
{
	auto TOT = cur_time();
	auto T = cur_time();

	eprn(":: Bucketing into {} files with {} threads", nfiles, nthreads);
	
	/// 1. load counts
	unordered_map<uint32_t, double> counts;
	string s, q, n, r;
	ifstream fin(path.c_str());
	while (getline(fin, s)) {
		uint32_t barcode = 0;
		DO(BC_LEN) barcode = (barcode << 2) | hash_dna(s[_]);
		assert(barcode != 0);
		counts[barcode] = 0;
	}
	fin.close();

	map<string, full_count_t> full_counts;
	glob_t glob_result;
	glob((prefix + "*.ema-count").c_str(), GLOB_TILDE, NULL, &glob_result);
	for (int fi = 0; fi < glob_result.gl_pathc; fi++) {
		string fp = glob_result.gl_pathv[fi];
		if (!S_ISREG(stat_file(fp))) 
			continue;
		FILE *f = fopen(fp.c_str(), "rb");
		eprnn(":: Loading counts file {} ... ", fp);
		int64_t total;
		fread(&total, 8, 1, f);
		fread(&total, 8, 1, f);
		eprnn(" {:n}/10x ", total);
		while (total--) {
			uint32_t bcd;
			fread(&bcd, 4, 1, f);
			int64_t cnt;
			fread(&cnt, 8, 1, f);
			counts[bcd] += cnt;
		}
		fread(&total, 8, 1, f);
		eprnn(" {:n}/total ", total);
		string b = string(BC_LEN, '#');
		while (total--) {
			fread(&b[0], 1, BC_LEN, f);
			int64_t cnt;
			fread(&cnt, 8, 1, f);
			full_counts[b].count += cnt;
		}
		fclose(f);
		eprn("done!");
	}
	double total_counts = 0;
	for (auto &c: counts) total_counts += c.second + 1;
	for (auto &c: counts) c.second = (c.second + 1) / total_counts;
	eprn(":: Loading counts ... done in {:.1f} s", elapsed(T)); T = cur_time();

	/// 2. process barcodes
	vector<thread> threads;
	vector<stats_t> stats(nthreads);
	for (int i = 0; i < nthreads; i++) {
		threads.push_back(thread(correct_barcode, 
			i, nthreads,
			ref(counts), ref(full_counts), ref(stats[i]), /*do_h2*/ false));
	}
	for (auto &t: threads) t.join(); 
	eprn(":: Fixing barcodes ... done in {:.1f} s", elapsed(T)); T = cur_time();
	for (int i = 1; i < nthreads; i++) {
		stats[0].total_passed += stats[i].total_passed;
		stats[0].total_passed_bcd += stats[i].total_passed_bcd;
		stats[0].not_bucketed += stats[i].not_bucketed;
		stats[0].not_bucketed_bcd += stats[i].not_bucketed_bcd;
		stats[0].corrected_h1 += stats[i].corrected_h1;
		stats[0].corrected_h1_bcd += stats[i].corrected_h1_bcd;
		stats[0].corrected_h2 += stats[i].corrected_h2;
		stats[0].corrected_h2_bcd += stats[i].corrected_h2_bcd;
	}
	eprn(":: Total processed: {:n} ({:n} barcodes)\n"
		 "        no barcode: {:n} ({:n} barcodes)\n"
		 "      H1-corrected: {:n} ({:n} barcodes)\n"
		 "      H2-corrected: {:n} ({:n} barcodes)", 
		 stats[0].total_passed, stats[0].total_passed_bcd, 
		 stats[0].not_bucketed, stats[0].not_bucketed_bcd,
		 stats[0].corrected_h1, stats[0].corrected_h1_bcd, 
		 stats[0].corrected_h2, stats[0].corrected_h2_bcd);
	// exit(0);

	/// 3. fix boundaries
	unordered_map<uint32_t, pair<int64_t, int>> file_offsets; // barcode -> (size, file)
	for (auto &b: full_counts) if (b.second.real_bcd != 0) {
		file_offsets[b.second.real_bcd].first += b.second.count;
	}

	vector<pair<int64_t, uint32_t>> fx;
	fx.reserve(file_offsets.size());
	for (auto &b: file_offsets) fx.push_back({b.second.first, b.first});
	sort(fx.begin(), fx.end(), greater<pair<int64_t, uint32_t>>());

	struct pf {
		int64_t size;
		int index;
		bool operator< (const pf &p) const { return size > p.size; }
	};
	priority_queue<pf> qq;
	vector<FILE*> fo;
	vector<pair<uint8_t*, size_t>> buffer;
	vector<int64_t> fo_cnt;
	const int BUFFER_SIZE = 10 * 1024 * 1024 + 10 * 1024;
	system(fmt::format("mkdir -p {}/ema", prefix).c_str());
	for (int i = 0; i < nfiles; i++) {
		qq.push(pf {0, i});
		FILE *f = fopen(fmt::format("{}/ema/bin-{:03}", prefix, i).c_str(), "wb");
		fo.push_back(f);
		fo_cnt.push_back(0);

		uint8_t *bf = (uint8_t*)malloc(BUFFER_SIZE); // 10MB + 10KB
		buffer.push_back({bf, 0}); // 10MB
	}
	{ 
		FILE *f = fopen(fmt::format("{}/ema/bin-nobc", prefix, nfiles).c_str(), "wb");
		fo.push_back(f);
		uint8_t *bf = (uint8_t*)malloc(BUFFER_SIZE); // 10MB + 10KB
		buffer.push_back({bf, 0}); // 10MB
	}
	for (auto &b: fx) {
		auto c = qq.top(); qq.pop();
		c.size += b.first;
		file_offsets[b.second].second = c.index;
		qq.push(c);
	}
	file_offsets[0] = {0, nfiles};
	eprn(":: File assignment ... done in {:.1f} s", elapsed(T)); T = cur_time();
	
	/// 4. assign stuff
	char bcd[BC_LEN + 1];
	bcd[BC_LEN] = 0;
	string b = string(BC_LEN, '#');
	eprnn(":: Writing {:n} barcodes ... ", file_offsets.size());
	while (getline(cin, n)) { 
		getline(cin, r);
		getline(cin, q);
		getline(cin, q);

		uint32_t barcode;
		DO(BC_LEN) {
			barcode = (barcode << 2) | hash_dna(s[_]);
			b[_] = hash_dna(r[_]) * QUAL_BASE + min(QUAL_BASE - 1, q[_] - ILLUMINA_QUAL_OFFSET);
		}
		
		if (counts.find(barcode) == counts.end()) {
			auto it = full_counts.find(b);
			if (it != full_counts.end()) {
				barcode = it->second.real_bcd;
			} else {
				barcode = 0;
			}
		}

		auto &ff = file_offsets[barcode];
		int file = ff.second;

		uint8_t *buff = buffer[file].first;
		size_t &buffi = buffer[file].second;

		if (barcode) {
			fo_cnt[file]++;
		}

		if (barcode != 0) {
			DO(BC_LEN) bcd[BC_LEN - _ - 1] = "ACGT"[barcode & 3], barcode >>= 2;
			memcpy(buff + buffi, bcd, BC_LEN); 
			buffi += BC_LEN;
			buff[buffi++] = ' ';
		}

		for (auto c: n) {
			if (isspace(c)) break; 
			buff[buffi++] = c;
		}
		buff[buffi++] = ' ';

		memcpy(buff + buffi, r.c_str() + BC_LEN + MATE1_TRIM, r.size() - BC_LEN - MATE1_TRIM); 
		buffi += r.size() - BC_LEN - MATE1_TRIM;	
			
		buff[buffi++] = ' ';

		memcpy(buff + buffi, q.c_str() + BC_LEN + MATE1_TRIM, q.size() - BC_LEN - MATE1_TRIM); 
		buffi += r.size() - BC_LEN - MATE1_TRIM;
		buff[buffi++] = ' ';
		
		getline(cin, s); 
		getline(cin, s); 
		memcpy(buff + buffi, s.c_str(), s.size()); 
		buffi += s.size();
		buff[buffi++] = ' ';
		getline(cin, s); 
		getline(cin, s); 
		memcpy(buff + buffi, s.c_str(), s.size()); 
		buffi += s.size();
		buff[buffi++] = '\n';

		if (buffi >= 10 * 1024 * 1024) {
			fwrite(buff, 1, buffi, fo[file]);
			buffi = 0;
		}
	}
	for (int fi = 0; fi < fo.size(); fi++) {
		fwrite(buffer[fi].first, 1, buffer[fi].second, fo[fi]);
		fclose(fo[fi]);
		free(buffer[fi].first);
	}
	eprn(" done in {:.1f} s", elapsed(T)); T = cur_time();
}

