double randomDouble(double min, double max) {	
	return (gsl_rng_uniform(gsl_ran_r) * (max - min)) + min;
}



//vararg function that does fscanf + checks if read ok etc.
//int readFile(



//Does what it says on the tin, file starts rewinded? :( can find current line in function, revert afterwards :)
int countLines(FILE *f) {
	rewind(f);
	int num_lines = 0;
	int ch = 0;	
	while(!feof(f)) {
		ch = fgetc(f);
		if(ch == '\n') {
			num_lines++;
		}
	}	
	rewind(f);
	return num_lines;
}



//to match fftw3 convention, using row-major format
//[x][y][z] = [z + cells_z * y + cells_z * cells_y * x]
//http://www.fftw.org/doc/Row_002dmajor-Format.html#Row_002dmajor-Format
//ONLY FOR 3D ARRAY OF SIZE cells_x * cells_y * cells_z
int arr_ind(int xx, int yy, int zz) {
	if (xx < 0 || xx >= cells[0] || yy < 0 || yy >= cells[1] || zz < 0 || zz >= cells[2]) printf("grid wtf \n");
	return zz + cells[2] * yy + cells[2] * cells[1] * xx;
}

//same as above, but for displacements array
int arr_ind_displ(int xx, int yy, int zz) {
	return zz + cells_displ[2] * yy + cells_displ[2] * cells_displ[1] * xx;
}

int NaN_check_double(double x, char err[]) {
	if (isnan(x) || x != x) {
		printf("NaN value. %s \n", err);
	}
	return 0;
}

char* prep_path(char tail[]) {
	char* toReturn = malloc(100 * sizeof(*toReturn));
	sprintf(toReturn, "%s%s", home_directory, tail);
	return toReturn;
}



double z_to_a(double z) {
	return 1.0/(z + 1.0); 
}

double a_to_z(double a) {
	return 1.0/a - 1.0;
}

double R_nonlinear(Spline* variance) {
	Spline *variance_reversed = reversed_spline(variance);
	double Rnl = splint_generic(variance_reversed, 1.0);
	free(variance_reversed);	
	return Rnl;
}

double H_a(double a, Parameters *params) {
	return params->H0*sqrt(params->omega_v_0 + params->omega_m_0*pow(a,-3.0) + params->omega_r_0*pow(a,-4.0) + (1.0 - params->omega)*pow(a,-2.0));	
}

double rho_bar_z(Parameters *params) {
	//rintf("rho cosm: %lf, z: %lf, a: %lf, omega m: %lf \n", rho_cosm, params->z, z_to_a(params->z), params->omega_m_0);
	return rho_cosm * pow(z_to_a(params->z), -3.0) * params->omega_m_0;
}

double m_to_R(double mass, Parameters *params) {
	return pow(3.0*mass/(4.0*pi*rho_bar_z(params)), 1.0/3.0);
}

double R_to_m(double R, Parameters *params) {
	//printf("R: %lf, rho bar: %lf, om m : %lf \n", R, rho_bar_z(params), params->omega_m_0);
	return (4.0/3.0)*pi*pow(R,3.0)*rho_bar_z(params);
}

double omega_m_z(double z, Parameters *params) {
	double a = z_to_a(z);
	double H = H_a(a, params);
	return params->omega_m_0*pow(a, -3.0)*pow(params->H0/H, 2.0);
}

//Bullock et al. 2006?? or 2004?
////https://arxiv.org/pdf/astro-ph/9908159v3.pdf
double del_nl_z(double z, Parameters *params) {
	double om_m = omega_m_z(z, params);
	return (18.0*pi*pi + 82.0*(om_m - 1.0) - 39.0*pow(om_m - 1.0, 2.0)) / om_m; 
}




/*
int print_anything(char outFile[], BinInfo *xbi, double (*y_x)(double)) {
	double x_cur, y_cur, dx;
	dx = (xmax - xmin) / (double)bins;
	FILE *f = fopen(outFile, "w");
	for (int i = 0; i <= bins; i++) {
		x_cur = xmin + (double)i *dx;
		y_cur = (*y_x)(x_cur);
		fprintf(f, "%le \t %le\n", x_cur, y_cur);	
	}
	fclose(f);
	return 0;
}
*/

void prep_out_path(char **holder, char folder[], char subfolder[], char filename[]) {
	char final_str[200];
	//sprintf(final_str, "%s%s%s

}

/*
int print_delsq(Spline *Pk_spline, char filename[], BinInfo* binInfo, bool scaled) {
	
	FILE *f = fopen(filename, "w");	
	double delsq, k;
		
	for (int i = 0; i < binInfo->bins; i++) {
		k = bin_to_x(binInfo, i);		
		if (scaled) {
			delsq = delsq_lin(Pk_spline, k*s);
		} else {
			delsq = delsq_lin(Pk_spline, k);
		}
		fprintf(f, "%lf \t %le", k, delsq);
		fprintf(f, "\n");
	}
	
	fclose(f);
	return 0;
}
*/
