#include "CubicSpline_double.c"
//
//
//MAIN FUNCTIONS
//
//

//Create new struct Spline with the spline coefficients, x values, y values, xmin, xmax, number of lines lines
//the xs must be monotonously increasing
Spline *new_spline(int lines, double *xs, double *ys) {
	
	double xmin = 0.0;
	double xmax = 0.0;	
	bool x_mono = true;
	bool x_inc = false;
	double *x_vals = malloc((size_t)lines * sizeof(*xs));
	double *y_vals = malloc((size_t)lines * sizeof(*ys));
	if ((x_vals != NULL) && (y_vals != NULL)) {
		memcpy(x_vals, xs, lines * sizeof(*x_vals));
		memcpy(y_vals, ys, lines * sizeof(*y_vals));
	} else {
		printf("malloc err in new_spline()\n");
		exit(0);
	}
	double *coeffs = malloc((size_t)lines * sizeof(*coeffs));	
	int dupcount = 0;
	int smalldiff_count = 0;		
	
	//monotoneity check for x
	for (int i = 0; i < lines-1; i++) {
		if (x_vals[i+1] > x_vals[i]) {
			x_inc = true;
		} else if (x_vals[i+1] < x_vals[i]) {
			//x decreased going forward, check monotonous
			if (x_inc == true) {
				//allow for error margin, if 'not monotonous'
				if ((x_vals[i] - x_vals[i+1])/x_vals[i] < 1e-9) {
					x_vals[i+1] = x_vals[i]*(1.0 + 1e-9); 
					smalldiff_count += 1;
				} else {
					printf("the x-values of the function given are not monotonous. %le followed by %le followed by %le \n", x_vals[i-1], x_vals[i], x_vals[i+1]);
					exit(0);
				}
			}
		} else {
			//x stayed the same -> nonsense??
			if (dupcount < 2) {
				printf("duplicate x-values in function given, check index %d. x[i]: %le, x[i+1]: %le\n", i, x_vals[i], x_vals[i+1]);
				dupcount += 1;	
			}			
		}
	}
	
	if (smalldiff_count > 0) printf("\nx values very close %d times\n\n", smalldiff_count);
		
	if (x_inc) {
		xmin = x_vals[0];
		xmax = x_vals[lines-1];
	} else {
		xmin = x_vals[lines-1];
		xmax = x_vals[0];
		//also reverse, spline/splint wants x in increasing order
		double temp, temp2;		
		for (int i = 0; i < lines/2; i++) {
			temp = x_vals[i];
			temp2 = y_vals[i];
			x_vals[i] = x_vals[lines - i - 1];
			y_vals[i] = y_vals[lines - i - 1];
			x_vals[lines - i - 1] = temp;
			y_vals[lines - i - 1] = temp2;
		}	
	}	
	
	//spline to get the coefficients, store spline information, together with coefficients, and return it	
	spline(x_vals, y_vals, lines, 1.0e31, 1.0e31, coeffs);		
	Spline *toReturn = malloc(sizeof(*toReturn));
	if (toReturn == NULL) perror("malloc toReturn failed in new_spline()\n");
	toReturn->x_vals = x_vals;
	toReturn->y_vals = y_vals;
	toReturn->coeffs = coeffs; 
	toReturn->lines = lines; 
	toReturn->xmin = xmin; 
	toReturn->xmax = xmax; 	
	return toReturn;	
}

Spline *copy_spline(Spline *s) {
	Spline *toReturn = malloc(sizeof(*toReturn));
	if (toReturn == NULL) perror("error in copy_spline(). ");

	toReturn->x_vals = malloc(s->lines * sizeof(*(s->x_vals)));
	if (toReturn->x_vals == NULL) {
		printf("not gud in copy_spline()\n");
	} else {
		memcpy(toReturn->x_vals, s->x_vals, s->lines*sizeof(*(s->x_vals)));
		if (!equalsDoubleArray(toReturn->x_vals, s->x_vals, s->lines)) {
			printf("bug in copy_spline()\n");
			exit(0);
		}
	}
	
	toReturn->y_vals = malloc(s->lines * sizeof(*(s->y_vals)));
	if (toReturn->y_vals == NULL) {
		printf("not gud in copy_spline()\n");
	} else {
		memcpy(toReturn->y_vals, s->y_vals, s->lines*sizeof(*(s->y_vals)));
		if (!equalsDoubleArray(toReturn->y_vals, s->y_vals, s->lines)) {
			printf("bug in copy_spline()\n");
			exit(0);
		}
	}
	
	toReturn->coeffs = malloc(s->lines * sizeof(*s->coeffs));
	if (toReturn->coeffs == NULL) {
		printf("not gud in copy_spline()\n");
	} else {
		memcpy(toReturn->coeffs, s->coeffs, s->lines*sizeof(*(s->coeffs)));
		if (!equalsDoubleArray(toReturn->coeffs, s->coeffs, s->lines)) {
			printf("bug in copy_spline()\n");
			exit(0);
		}
	}
	
	toReturn->lines = s->lines;
	toReturn->xmin = s->xmin;
	toReturn->xmax = s->xmax;	
	return toReturn;
}

int free_spline(Spline *s) {
	free(s->x_vals);
	free(s->y_vals);
	free(s->coeffs);
	free(s);
	s = NULL;
	return 0;
}

//allowing small deviation to counter bugs with the result ending up 0 inside the range
double splint_generic(Spline *spline, double x) {
	double y = 0.0;
	double margin = (spline->xmax - spline->xmin)*1e-6;	
	//printf("splinting generic, x: %lf\n", x);
	//print_spline(spline, NULL);
	if ((x > spline->xmin - margin) && (x < spline->xmax + margin)) {
		//printf("splint in range\n");		
		if (spline->xmin - margin < x && x < spline->xmin + margin) {
			y = spline->y_vals[0];
		} else if (spline->xmax - margin < x && x < spline->xmax + margin) {
			y = spline->y_vals[spline->lines - 1];
		} else {
			splint(spline->x_vals, spline->y_vals, spline->coeffs, spline->lines, x, &y);
		}	
	}	
	return y;
}

int print_spline(Spline *spline, BinInfo *binInfo) {
	printf("PRINTING SPLINE------------------------\n");
	//just printing the arrays
	if (binInfo == NULL) {
		for (int i = 0; i < spline->lines; i++) {
			printf("x: %le, y: %le\n", spline->x_vals[i], spline->y_vals[i]);
		}
	} else {
		double x = 0.0;
		for (int i = 0; i < binInfo->bins; i++) {			
			x = bin_to_x(binInfo, i);			
			printf("x: %le, y: %le\n", x, splint_generic(spline, x)); 
		}		
	}
	printf("DONE PRINTING SPLINE-------------------\n");
	return 0;
}

int print_spline_file(Spline *spline, BinInfo *binInfo, char *outFile) {
	FILE *f = fopen(outFile, "w");
	for (int i = 0; i < binInfo->bins; i++) {
		fprintf(f, "%le \t %le", bin_to_x(binInfo, i), splint_generic(spline, bin_to_x(binInfo, i))); 
		fprintf(f, "\n");
	}
	fclose(f);
	return 0;

}


//
//
//OTHER FUNCTIONS
//
//

Spline *input_spline_file(char *filename, char *format) {
	FILE *f = fopen(filename, "r");
	int lines = countLines(f);
	double *xs = malloc((size_t)lines * sizeof(*xs));
	double *ys = malloc((size_t)lines * sizeof(*ys));
	if ((xs != NULL) && (ys != NULL)) {
		int par = 0;
		for (int i = 0; i < lines; i++) {	
			par = fscanf(f, format, &xs[i], &ys[i]);
			if (par != 2) {
				printf("wrong spline file format? input_spline_file()\n");
				exit(0);
			}	
		}
		Spline *toReturn = new_spline(lines, xs, ys);
		free(xs);
		free(ys);
		xs = NULL;
		ys = NULL;
		printf("input spline from %s with %d lines, xmin: %le, xmax: %le\n", filename, toReturn->lines, toReturn->xmin, toReturn->xmax);	
		return toReturn;
	} else {
		printf("error in input_spline_file()\n");
		exit(0);	
		return NULL;	
	}	
}
/*
//changes f(x) -> af(x)
//should actually times coeffs by scale_fac? or not that simple?
int scale_spline(Spline **spline, double scale_fac) {
	double *ys = malloc((*spline)->lines * sizeof(*ys));
	if (ys == NULL) {
		printf("err in scale_spline()\n");
		exit(0);
	}
	for (int i = 0; i < (*spline)->lines; i++) {
		ys[i] = scale_fac * (*spline)->y_vals[i];
	}
	Spline *toReturn = new_spline((*spline)->lines, (*spline)->x_vals, ys);
	free_spline(*spline);
	free(ys);
	*spline = toReturn;
	return 0;
}*/

//changes f(x) -> af(x)
//should actually times coeffs by scale_fac? or not that simple?
int scale_spline(Spline *s, double scale_fac) {
	double *ys = malloc(s->lines * sizeof(*ys));
	double *coeffs = malloc(s->lines * sizeof(*coeffs));
	if (ys == NULL) {
		printf("err in scale_spline()\n");
		exit(0);
	}
	for (int i = 0; i < s->lines; i++) {
		ys[i] = scale_fac * s->y_vals[i];
	}
	spline(s->x_vals, ys, s->lines, 1.0e31, 1.0e31, coeffs);
	free(s->y_vals);
	s->y_vals = ys;
	free(s->coeffs);
	s->coeffs = coeffs;	
	return 0;
}

//FIX AS ABOVE
int spline_add_const(Spline *s, double c) {
	double *ys = malloc(s->lines * sizeof(*ys));
	double *coeffs = malloc(s->lines * sizeof(*coeffs));
	for (int i = 0; i < s->lines; i++) {
		ys[i] = s->y_vals[i] + c;
	}
	spline(s->x_vals, ys, s->lines, 1.0e31, 1.0e31, coeffs);
	free(s->y_vals);
	s->y_vals = ys;
	free(s->coeffs);
	s->coeffs = coeffs;	
	return 0;
}


/*
//f(x) -> f(ax)
Spline *scale_x(Spline *spline, double sf) {
	double *xs = malloc(spline->splineInfo->lines * sizeof(*xs));
	for (int i = 0; i < spline->splineInfo->lines; i++) {
		xs[i] = spline->splineInfo->x_vals[i] / sf;		
	}
	Spline *toReturn = input_spline_values(spline->splineInfo->lines, xs, spline->splineInfo->y_vals);
	free(spline);
	free(xs);	
	return toReturn;
}


*/


//sets the value at the pointer "reversed" to the reverse x = f(y) function of the spline "current"
Spline *reversed_spline(Spline *current) {
	Spline *reversed = new_spline(current->lines, current->y_vals, current->x_vals);
	return reversed;
}


//prepares a spline from a function of the type double y(double x){}
Spline *prep_spline_generic(BinInfo *binInfo, int reverse, double (*func)(double)) {
	double x_cur = 0.0;
	double y_cur = 0.0;
	Spline *toReturn = NULL;
	double *xs = malloc((size_t)binInfo->bins * sizeof(*xs));
	double *ys = malloc((size_t)binInfo->bins * sizeof(*ys));	
	
	for (int i = 0; i < binInfo->bins; i++) {
		x_cur = bin_to_x(binInfo, i);
		y_cur = (*func)(x_cur);
		xs[i] = x_cur;
		ys[i] = y_cur;
	}	
	
	if (reverse == REVERSED) {
		toReturn = new_spline(binInfo->bins, ys, xs);
	} else if (reverse == NORMAL) {
		toReturn = new_spline(binInfo->bins, xs, ys);
	} else {
		printf("prep_spline_generic expected REVERSED or NORMAL \n");
		exit(0);
	}	
	return toReturn;
}







