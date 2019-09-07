//////////////////////////////////////////////////////////////////////
// main.cpp -- Test Program for Plan13 class
///////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "plan13.hpp"

#include <array>

static std::array<const char*[3],3> elements({{
	{"ISS",
	 "1 25544U 98067A   19247.71039296  .00002175  00000-0  45453-4 0  9990",
	 "2 25544  51.6465 326.8109 0008020   5.4414 144.2903 15.50432168187594"},
	{"AO-7",
	 "1 07530U 74089B   19247.19524406 -.00000034  00000-0  67221-4 0  9997",
	 "2 07530 101.7517 214.0840 0011902 196.0911 234.8002 12.53639545 50230"},
	{"UO-11",
	 "1 14781U 84021B   19247.49076266  .00000042  00000-0  11238-4 0  9994",
	 "2 14781  97.6079 268.5149 0009924  68.5131 291.7140 14.83103747910297"},
}});

static void
write_func(const void *buf,unsigned bytes,void *arg) {
	(void)arg;	// Not used here

	::write(1,buf,size_t(bytes));	// Write to stdout
}

int
main(int argc,char **argv) {
	Plan13 p13(write_func);

	// p13.set_frequency(435300000,1459200000);
	p13.set_location(-79.000, 43.1100, 100); // QTH
	p13.init();

	p13.print_header();
	for (unsigned x = 0; x < elements.size(); x++) {
		p13.load_elements(elements[x][0],elements[x][1],elements[x][2]);
		p13.set_time(time(nullptr));
		p13.calculate(); 
		p13.print();
	}
}

// End main.cpp
