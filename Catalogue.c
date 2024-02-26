int free_catalogue(Catalogue *cat) {
	for (int i = 0; i < cat->particle_no; i++) {
		if (cat->particles[i] != NULL) {
			free_particle(cat->particles[i]);
		} else {
			printf("trying to free catalogue with particles NULL inside? \n");
			exit(0);
		}
	}
	free(cat->particles);
	cat->particles = NULL;
	free(cat);
	cat = NULL;
	return 0;
}

Catalogue *copy_catalogue(Catalogue *src) {
	Catalogue *toReturn = malloc(sizeof(*toReturn));

	Particle **particles = malloc((size_t)src->particle_no * sizeof(*particles));
	for (int i = 0; i < src->particle_no; i++) {
		particles[i] = copy_particle(src->particles[i]);
	}
	
	toReturn->particle_no = src->particle_no;
	toReturn->particles = particles;
	return toReturn;
}

Catalogue *combine_catalogues(Catalogue *cat1, Catalogue *cat2) {
	Catalogue *toReturn = malloc(sizeof(*toReturn));
	int total_particles = cat1->particle_no + cat2->particle_no;
	Particle **particles = malloc((size_t)total_particles * sizeof(*particles));
	
	//deep copy
	for (int i = 0; i < cat1->particle_no; i++) {		
		particles[i] = copy_particle(cat1->particles[i]);		
	}
	
	for (int i = 0; i < cat2->particle_no; i++) {
		particles[i + cat1->particle_no] = copy_particle(cat2->particles[i]);
	}
	
	printf("\nCombined catalogues with %d and %d particles, final catalogue %d particles\n", cat1->particle_no, cat2->particle_no, total_particles);

	toReturn->particles = particles;
	toReturn->particle_no = total_particles;
	return toReturn;
}

int catalogue_to_file(Catalogue *cat, char outFile[]) {
	FILE *f = fopen(outFile, "w");
	if (f == NULL) {
		printf("no such file %s \n", outFile);
		exit(0);
	} 
	for (int i = 0; i < cat->particle_no; i++) {
		fprintf(f, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %le %d\n", cat->particles[i]->x, cat->particles[i]->y, cat->particles[i]->z, cat->particles[i]->vx, cat->particles[i]->vy, cat->particles[i]->vz, cat->particles[i]->mass, cat->particles[i]->parent_id);
	}
	fclose(f);
	return 0;
}

//finds the min and max mass of a catalogue
double *mass_min_max(Catalogue *catalogue) {
	double Mmin, Mmax;
	Mmax = Mmin = catalogue->particles[0]->mass;
	double *toReturn = malloc(2*sizeof(*toReturn));	
	for (int i = 1; i < catalogue->particle_no; i++) {
		if (catalogue->particles[i]->mass < Mmin) {
			Mmin = catalogue->particles[i]->mass;
		} else if (catalogue->particles[i]->mass > Mmax) {
			Mmax = catalogue->particles[i]->mass;
		}
	}
	toReturn[0] = Mmin;
	toReturn[1] = Mmax;
	return toReturn;
}

int catalogue_setMasses(Catalogue *cat, Spline *cumsums_reverse) {
	double cumsum_min = cumsums_reverse->xmin;
	double cumsum_max = cumsums_reverse->xmax;
	double rand_cumsum, cur_mass;	
	for (int i = 0; i < cat->particle_no; i++) {
		rand_cumsum = randomDouble(cumsum_min, cumsum_max);
		cur_mass = splint_generic(cumsums_reverse, rand_cumsum);
		if (cur_mass < 1e-6 || cur_mass < splint_generic(cumsums_reverse, cumsum_min)) {
			printf("cur mass: %le, rand_cumsum: %le, cumsum min: %le, cumsum max: %le in catalogue_setMasses. ZA mocks?\n", cur_mass, rand_cumsum, cumsum_min, cumsum_max);
			exit(0);
		}
		cat->particles[i]->mass = cur_mass;
	}
	double *mass_lims = mass_min_max(cat);
	printf("set masses according to cumulative sum reverse spline given. Mmin: %le, Mmax: %le\n", mass_lims[0], mass_lims[1]);
	free(mass_lims);
	return 0;
}

int randomise_catalogue(Catalogue *cat) {	
	int ind1 = 0;
	int ind2 = 0;		
	for (int i = 0; i < 2*cat->particle_no; i++) {
		ind1 = (int)randomDouble(0.0, (double)cat->particle_no);
		ind2 = (int)randomDouble(0.0, (double)cat->particle_no);
		printf("rnd %d\n", i);
		Particle *tmp = copy_particle(cat->particles[ind1]);
		free_particle(cat->particles[ind1]);
		cat->particles[ind1] = copy_particle(cat->particles[ind2]);
		free_particle(cat->particles[ind2]);
		cat->particles[ind2] = tmp;
	}
	return 0;
}

Catalogue *input_catalogue_file(const char filename[], int skipLines, const char format[]) {

	FILE *catalogue_input = fopen(filename, "r");	
	Catalogue *this_cat;
	int fs = 0;
	
	if (catalogue_input == NULL) {
		printf("\nerror opening %s in input_catalogue_file()\n", filename);
		exit(0);
	} else {
		//count lines
		int num_lines = countLines(catalogue_input);
	
		//ALLOCATE MEMORY TO PARTICLES
		int particle_no = num_lines - skipLines;
		Particle **particles = malloc((size_t)particle_no * sizeof(*particles));	
		
		//skip first [skipLines] lines, due to the format of file
		for (int i = 0; i < skipLines; i++) {
			fs = fscanf(catalogue_input, "%*[^\n]\n");
		}
		
		//create catalogue
		this_cat = malloc(sizeof(*this_cat));
		
		//x, y, z in (Mpc/h), then vx, vy, vz in (km/s), then mass in (solar masses); \\\\\\%* ignore flag		
		for (int i = 0; i < particle_no; i++) {
		
			Particle *particle = malloc(sizeof(*particle));
			fs = fscanf(catalogue_input, format, &particle->x, &particle->y, &particle->z,
			&particle->vx, &particle->vy, &particle->vz, &particle->mass);			
			
			char *label	= malloc(20 * sizeof(*label));
			sprintf(label, "particle %d", i);
			particle->label = label;
			particle->type = UNKNOWN;		
			
			if (fs != 7) {
				printf("format in input catalogue mismatch? Read %d values in row %d. Format: %s\n", ph, i, format);
				exit(0);
			}
			particles[i] = particle;
			//next line
			fs = fscanf(catalogue_input, "%*[^\n]\n");			
		}
		
		this_cat->particle_no = particle_no;
		this_cat->particles = particles;
		
		printf("\ncatalogue from %s input successfully. Particles: %d, using a box of size %lf (IS THIS THE RIGHT SIZE?) \n", filename, particle_no, volume_limits[0]);	
	}
	
	fclose(catalogue_input);
	return this_cat;
}

