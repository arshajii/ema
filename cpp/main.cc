/// 786

#include <iostream>
#include <string>

#include "common.h"
#include "count.h"
#include "correct.h"

using namespace std;

/******************************************************************************/

int main(int argc, char **argv) 
{
	ios_base::sync_with_stdio(0);
	setlocale(LC_NUMERIC, "en_US.UTF-8");

	string bcd_path = "4M-with-alts-february-2016.txt";

	string cmd = argv[1];
	if (cmd == "count") {
		count(bcd_path);
	} else if (cmd == "correct") {
		// prefix nthreads nfiles
		correct(bcd_path, argv[2], atoi(argv[3]), atoi(argv[4]));
	} else {
		eprn("whoops");
	}

	return 0;
}
