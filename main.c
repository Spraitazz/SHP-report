#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
//#include <mcheck.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <fftw3.h>

//my definitions
#define pi 3.14159265358979323846264338327

#define X 0
#define Y 1
#define Z 2

#define REAL_SPACE 2221
#define REDSHIFT_SPACE 2222

#define NGP 3331
#define CIC 3332

#define NORMAL 4441
#define REVERSED 4442

#define NORMAL_BIN 5551
#define LOG_BIN 5552

#define F4 1111
#define F5 1112
#define F6 1113
#define GR 1114

#define MEASURED 6661
#define MODEL 6662

#define HALO 7771
#define SATELLITE 7772
#define CENTRAL 7773
#define UNKNOWN 7774

//my headers
#include "structs.h"
#include "main.h"

#include "basic_func.c"

//"objects"
#include "BinInfo.c"
#include "Spline.c"
#include "Pk_Spline.c"
#include "Particle.c"
#include "Catalogue.c"

//my code
#include "memory.c"
#include "parameters.c"
#include "functions.c"
#include "power_spectrum.c"
//#include "fold.c"
#include "runs.c"
//#include "HOD.c"
//#include "HOD_tests.c" //should be running HOD.c
//#include "covariance_matrices.c"
#include "bias.c"
#include "ZA.c"
#include "Dplus.c"
#include "rescaling_functions.c"
#include "mocks.c"
//#include "rescaling.c"
//#include "RSD_Pk.c"
//#include "KaiserGauss.c"
//#include "halo_model.c"
//#include "Pk_model.c"

#include "make_mocks.c"

//#include "tests.c"

//http://arxiv.org/pdf/astro-ph/0005010v2.pdf rescaling
//http://arxiv.org/pdf/1408.1047v2.pdf rescaling
//http://arxiv.org/pdf/1308.5183v3.pdf rescaling
//http://arxiv.org/pdf/1202.5559v4.pdf HOD
//http://arxiv.org/pdf/1005.2413v2.pdf HOD parameters


//http://halotools.readthedocs.io


int main(int argc, char **argv) {	

	//mcheck(MCHECK_HEAD);
	
	//unused parameters
	argc = argc;
	argv = argv;	
	//run time
	clock_t start_time = clock();

	//set home directory!!
	//home_directory = "/shome/jonasp/Work/Testing_GR/code/Jonas";	
	home_directory = "/home/jonas/Testing_GR/code/Jonas";
	homeDirChk();
	sprintf(temp_directory, "%s/output/tmp", home_directory);	
	
	//define volume for particles
	volume_limits[0] = volume_limits[1] = volume_limits[2] = 256.0;  // h^-1 Mpc 	
	volume = volume_limits[0]*volume_limits[1]*volume_limits[2];	
	
	//define P(k) FFT grid size, EVEN NUMBERS ONLY TO AVOID PROBLEMS, 2^n best for FFT	
	cells[0] = cells[1] = cells[2] = 128;	
	
	printf("Expect folding to start around sqrt(3)*k_nyq/4 = %lf \n", (sqrt(3.0)*pi / ((volume_limits[0] / (double) cells[0]))/4.0));
	
	//grid size for ZA displacements
	cells_displ[0] = 128;
	cells_displ[1] = 128;
	cells_displ[2] = 128;	
		
	//spectrum bins
	spectrum_size = 30;	
		
	//start random numb generator
	init_rng();			
	
	//HOD RANDS STUFF	
	set_params();
	
	
	//covariance_sequence();
	//HOD_main();
	//biased_displacements();
	//ZA_RSD();
	//ST_bias();
	
	//HOD_mocks_test2();
	//final_rescaling();
	
	make_mocks();
	measure_Pk();
	
	//rescale_randoms_model_to_model_redshiftspace();
	//HOD_test();
	//ZA_mocks_HOD_test();
	
	//print_cdm("/home/jonas/Testing_GR/code/Jonas/data/cdm.dat");

	//rescale_testRun();
	//rescale_model();
	//rescale_catalogue_to_model();
	//rescale_catalogue_to_catalogue();
	exit(0);
	
	
	
	
	clock_t end_time = clock();
	printf("\nDone in %.1lf seconds\n", (double)(end_time - start_time)/CLOCKS_PER_SEC);

	return 0;
	
	
	
	
	
}
