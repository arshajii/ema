/// 786

#include <iostream>
#include <string>

#include "common.h"
#include "count.h"
#include "correct.h"

using namespace std;

/******************************************************************************/

EXTERNC void cppinit()
{
	ios_base::sync_with_stdio(0);
	setlocale(LC_NUMERIC, "en_US.UTF-8");
}

int mainx(int argc, const char **argv)
{
	cppinit();

	const char *bcd_path = "4M-with-alts-february-2016.txt";

	string cmd = argv[1];
	if (cmd == "count") {
		count(bcd_path, argv[2], 1 * GB);
	} else if (cmd == "correct") {
		// prefix nthreads nfiles
		const char *output_dir = argv[2];
		int nthreads = atoi(argv[3]);
		int nbuckets = atoi(argv[4]);
		correct(bcd_path,
			argv + 5,
			argc - 5,
			output_dir,
			0,
			10 * MB,
			nthreads,
			nbuckets);
	} else {
		eprn("whoops");
	}

	return 0;
}
