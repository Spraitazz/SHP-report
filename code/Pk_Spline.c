//Why does P(k) care about H0? Maybe just have redshift??
//make Spline the bare class, Pk_Spline on top - has model, has redshift
Pk_Spline *make_Pk_spline(Spline *spline, double z) {
	Pk_Spline *toReturn = malloc(sizeof(*toReturn));
	if (toReturn == NULL) perror("error in make_Pk_spline(). ");	
	toReturn->spline = copy_spline(spline);
	toReturn->pk_loA = 0.0;
	toReturn->pk_hiA = 0.0;
	toReturn->pk_lon = 0.0;
	toReturn->pk_hin = 0.0;
	toReturn->model = MEASURED;
	toReturn->z = z;
	return toReturn;
}

Pk_Spline *copy_Pk_spline(Pk_Spline *pk_s) {
	Pk_Spline *pk_cpy = make_Pk_spline(pk_s->spline, pk_s->z);
	//print_spline(pk_cpy->spline, NULL);
	pk_cpy->pk_loA = pk_s->pk_loA;
	pk_cpy->pk_hiA = pk_s->pk_hiA;
	pk_cpy->pk_lon = pk_s->pk_lon;
	pk_cpy->pk_hin = pk_s->pk_hin;
	pk_cpy->model = pk_s->model;
	pk_cpy->z = pk_s->z;	
	return pk_cpy;
}

int free_Pk_spline(Pk_Spline *pk_spline) {	
	free_spline(pk_spline->spline);
	free(pk_spline);
	pk_spline = NULL;
	return 0;
}


//add power low regression tails to both ends of spline
int extend_spline(Pk_Spline *pk_spline) {	

	double pk_loA = 0.0;
	double pk_hiA = 0.0;
	double pk_lon = 0.0;
	double pk_hin = 0.0;
	double xmin = pk_spline->spline->xmin;
	double xmax = pk_spline->spline->xmax;	
	double xend_min = exp(log(xmin) + 0.2 * (log(xmax) - log(xmin))); //10% of range, can vary this
	double xstart_max = xmax - 0.1 * (xmax - xmin);	
	
	powerlaw_regression(pk_spline->spline->lines, xmin, xend_min, 1.0, pk_spline->spline->x_vals, pk_spline->spline->y_vals, &pk_loA, &pk_lon);    
    powerlaw_regression(pk_spline->spline->lines, xstart_max, xmax, 1.0, pk_spline->spline->x_vals, pk_spline->spline->y_vals, &pk_hiA, &pk_hin);    

	pk_spline->pk_hiA = pk_hiA;
	pk_spline->pk_loA = pk_loA;
	pk_spline->pk_hin = pk_hin;
	pk_spline->pk_lon = pk_lon;
	pk_spline->model = MODEL;
	return 0;
}


//uses spline from what's available + power laws at low and high k. (all here return VP(k))
double splint_Pk_model(Pk_Spline *pk_model, double k) {
	double Pk = 0.0;
	if (pk_model->model != MODEL) {
		printf("asking to splint a model, but giving a non-model extended SplineInfo\n");
		exit(0);
	} else {
		if (k < pk_model->spline->xmin) {
			Pk = pk_model->pk_loA*pow(k, 3.0 + pk_model->pk_lon);		
		} else if (k > pk_model->spline->xmax) {
			Pk = pk_model->pk_hiA*pow(k, 3.0 + pk_model->pk_hin);		
		} else {
			Pk = splint_generic(pk_model->spline, k);
		}
	}
	return Pk/volume;
}

double splint_Pk(Pk_Spline *pk_spline, double k) {
	if (pk_spline->model == MODEL) {
		return splint_Pk_model(pk_spline, k);
	} else {
		//printf("splinting nonmodel\n");
		return splint_generic(pk_spline->spline, k)/volume;
	}
}

Pk_Spline *input_Pk_spline_file(char *inFile, char *format, double z) {
	Spline *s = input_spline_file(inFile, format);
	Pk_Spline *pk_s = make_Pk_spline(s, z);
	free_spline(s);
	return pk_s;
}

int print_Pk_spline(Pk_Spline *pk_spline, BinInfo *binInfo) {
	print_spline(pk_spline->spline, binInfo);
	return 0;
}

int print_Pk_spline_file(Pk_Spline *pk_s, BinInfo *binInfo, char *outFile) {
	print_spline_file(pk_s->spline, binInfo, outFile);
	return 0;
}
