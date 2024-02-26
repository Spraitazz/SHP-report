Catalogue *Jenkins_fold_volume(Catalogue *catalogue) {	

	//copy before changing
	Catalogue *toReturn = copy_catalogue(catalogue);

	volume_limits[0] /= Jenkins_foldfactor;
	volume_limits[1] /= Jenkins_foldfactor;
	volume_limits[2] /= Jenkins_foldfactor;			
	
	for (int i = 0; i < catalogue->particle_no; i++) {	
		catalogue->particles[i]->x = fmod(catalogue->particles[i]->x, volume_limits[0]);
		catalogue->particles[i]->y = fmod(catalogue->particles[i]->y, volume_limits[1]);
		catalogue->particles[i]->z = fmod(catalogue->particles[i]->z, volume_limits[2]);
	}
	return toReturn;
}

int Jenkins_restore(Catalogue *catalogue, Catalogue *catalogue_saved) {
	
	//free_catalogue(catalogue);
	for (int i = 0; i < catalogue->particle_no; i++) {
		catalogue->particles[i]->x = catalogue_saved->particles[i]->x;
		catalogue->particles[i]->y = catalogue_saved->particles[i]->y;
		catalogue->particles[i]->z = catalogue_saved->particles[i]->z;
	}
	//catalogue = copy_catalogue(catalogue_saved);
	free_catalogue(catalogue_saved);	
	
	volume_limits[0] *= Jenkins_foldfactor;
	volume_limits[1] *= Jenkins_foldfactor;
	volume_limits[2] *= Jenkins_foldfactor;		
	Jenkins_foldfactor = 1.0;	

	return 0;
}

FoldingInformation* combine_folded(const char unfolded_path[], const char folded_path[], char pk_format[], char outFile[], double kmin, double kmax, double k_fold, int* foldPoints, bool print) {

	FILE *unfolded_file = fopen(unfolded_path, "r");
	FILE *folded_file = fopen(folded_path, "r");
	FILE *output = fopen(outFile, "w");	
	
	double min_diff, k_diff, k_unfolded, k_folded, mono_unfolded;
	double mono_folded, perc_diff_k, perc_diff_mono;
	
	int len_folded, len_unfolded, min_diff_ind;
	int uf_first_ind, uf_last_ind, uf_alias_ind; 
	int f_first_ind, f_last_ind;
	int v = 0;	
	
	bool max_reached = false;
	bool success = true;
	
	FoldingInformation *toReturn = malloc(sizeof(*toReturn));

	//CHECK FOR NULL FILES
	if (unfolded_file == NULL) {
		printf("error opening %s\n", unfolded_path); 
		exit(0);
	}
	if (folded_file == NULL) {
		printf("error opening %s\n", folded_path);
		exit(0);
	}
	if (output == NULL) {
		printf("error opening %s\n", outFile);
		exit(0);
	}
	if (kmin < 0.0) {	
		printf("kmin can only be more or equal to 0.0\n");
		exit(0);	
	}
	
	//count lines
	len_folded = countLines(folded_file);
	len_unfolded = countLines(unfolded_file);
	//printf("len uf: %d, len f: %d \n", len_unfolded, len_folded);

	//[k][mono][quadro]
	double** folded_values;
	double** unfolded_values;
	calloc2D_double(&folded_values, len_folded, 6);
	calloc2D_double(&unfolded_values, len_unfolded, 6);	
	int folded_kBins[len_folded];
	int unfolded_kBins[len_unfolded];
	
	//read into arrays from the respective files
	for (int i = 0; i < len_folded; i++) {
		ph = fscanf(folded_file, pk_format, &folded_values[i][0], &folded_values[i][1], &folded_values[i][2], &folded_values[i][3], &folded_kBins[i], &folded_values[i][4], &folded_values[i][5]);		
	}	
	
	for (int i = 0; i < len_unfolded; i++) {
		ph = fscanf(unfolded_file, pk_format, &unfolded_values[i][0], &unfolded_values[i][1], &unfolded_values[i][2], &unfolded_values[i][3], &unfolded_kBins[i], &unfolded_values[i][4], &unfolded_values[i][5]);
	}	
		
/*
	if (kmax > folded_values[len_folded - 1][0]) {
		printf("kmax: %lf, more than maximum in folded file, %lf \n", kmax, folded_values[len_folded - 1][0]);	
		//exit(0);	
	}
	*/
	
	
	//-----------------------------------------------------------------------------------------
	
	//IMPORTANT: indices for connecting unfolded to folded		

	if (foldPoints == NULL) {
		//printf("foldpoints null \n");
	
		//find the k index from which to begin in the unfolded file	
		v = 0;
		while (unfolded_values[v][0] < kmin) {
			v+=1;
		}
		uf_first_ind = v;	
		
		//unfolded file start of aliasing
		v = len_unfolded - 1;	
		while (unfolded_values[v-1][1] < unfolded_values[v][1]) {
			v -= 1;
		}
		uf_alias_ind = v;
		
		//folded file "start of aliasing"
		v = 0;
		v = len_folded - 1;	
		while (folded_values[v-1][1] < folded_values[v][1] || folded_values[v-2][1] < folded_values[v-1][1]) {
			//printf("v: %d \n", v);
			v -= 1;
		}
		f_last_ind = v;	

		//maximum is input, check for finish of file before aliasing kicks in
		if (kmax > 0.0) {
			v = 0;
			while (folded_values[v][0] < kmax) {
				v += 1;
			}

			//finish before aliasing
			if (v < f_last_ind) {			
				f_last_ind = v;
				max_reached = true;
			}
		}

	
		//--------------------------------------------------------
	
		if (k_fold < 0.0) {
			//finding fold point automatically		
	
			//find the minimum k and monopole differences combining pairs from the folded/unfolded files	
			min_diff = DBL_MAX;

			//looping through unfolded file
			for (int i = uf_first_ind; i < uf_alias_ind; i++) {
				k_unfolded = unfolded_values[i][0];
				mono_unfolded = unfolded_values[i][1];
			
				//folded file
				for (int j = 0; j < f_last_ind; j++) {
					k_folded = folded_values[j][0];
					mono_folded = folded_values[j][1];
				
					//make sure k is going up from unfolded to folded
					if (k_folded > k_unfolded) {
						perc_diff_k = (k_folded - k_unfolded) / k_unfolded;
						perc_diff_mono = fabs((mono_folded - mono_unfolded) / mono_unfolded);
					
						//could POTENTIALLY just make sure monopole decreases with increasing k here (now more general)
						if (perc_diff_mono*perc_diff_k < min_diff && perc_diff_k < 0.2) {
											
							min_diff = perc_diff_mono*perc_diff_k;
							//best connection indices
							uf_last_ind = i;
							f_first_ind = j;						
						
						} 
					}
				}				
			}	
		} else {
		
			//folding point set, just find the necessary indices for connecting
			if (k_fold < kmin || k_fold > kmax) {
				printf("cannot fold below kmin OR above kmax\n");
				exit(0);
			} else {
				//just find the last index where k < k_fold in the unfolded file,
				//then jump to where k > k_fold in the folded file
				v = 0;
				while (unfolded_values[v][0] < k_fold) {
					v+=1;
				}
				uf_last_ind = v;
			
				v = 0;
				while(folded_values[v][0] < k_fold) {
					v+=1;
				}
				f_first_ind = v;		
			}	
		}
	
	} else {
	
		uf_first_ind = foldPoints[0];
		uf_last_ind = foldPoints[1];
		f_first_ind = foldPoints[2];
		f_last_ind = foldPoints[3];
		
		v = 0;
		while (folded_values[v][0] < kmax) {
			v += 1;
		}
	
		//finish before aliasing
		if (v < f_last_ind) {			
			f_last_ind = v;
			max_reached = true;
		}
		
		if (uf_last_ind >= len_unfolded || f_last_ind >= len_folded) {
			printf("file length mismatch \n");
			fclose(unfolded_file);
			fclose(folded_file);
			fclose(output);
			free2D_double(&folded_values, len_folded);
			free2D_double(&unfolded_values, len_unfolded);
			toReturn->success = false;
			return toReturn;						
		}
		
		//printf("in folding(), uf first ind: %d, uf last ind: %d, f first ind: %d, f last ind: %d \n", uf_first_ind, uf_last_ind, f_first_ind, f_last_ind);
		
	}	
	
	int lines_out = (uf_last_ind - uf_first_ind) + (f_last_ind - f_first_ind);
	
	//INFO TO SEND TO MULTI-FILE FUNCTION
	/*
	FoldingInformation toReturn = {
	.folded_values = folded_values,
	.unfolded_values = unfolded_values,
	.kmin = unfolded_values[min_k_ind_unfolded][0],
	.kmax = folded_values[aliasing_index_folded-1][0],
	.unfold_last_ind = unfold_last_ind,
	.fold_first_ind = fold_conn_ind,
	.min_k_ind = min_k_ind_unfolded,
	.lines = lines_out,
	.aliasing_index = aliasing_index_folded,
	.max_reached = max_reached
	};
	*/
	
	
	toReturn->uf_first_ind = uf_first_ind;
	toReturn->uf_last_ind = uf_last_ind;
	toReturn->f_first_ind = f_first_ind;
	toReturn->f_last_ind = f_last_ind;
	toReturn->max_reached = max_reached; 
	toReturn->success = true;
	
	//-----------------------------------------------------------------------	

	//if printing a single example
	if (print) {
		
		double final_values[lines_out][6];
		int final_kBins[lines_out];
	
		for (int i = uf_first_ind; i < uf_last_ind; i++) {
			for (int k = 0; k < 6; k++) {
				final_values[i - uf_first_ind][k] = unfolded_values[i][k];				
			}
			final_kBins[i - uf_first_ind] = unfolded_kBins[i];
		}
		//printf("done unfolded \n");		

		for (int i = f_first_ind; i < f_last_ind; i++) {					
			for (int k = 0; k < 6; k++) {
				final_values[(uf_last_ind - uf_first_ind) + i - f_first_ind][k] = folded_values[i][k];	
			}
			final_kBins[(uf_last_ind - uf_first_ind) + i - f_first_ind] = folded_kBins[i];
		}

		for (int i = 0; i < lines_out; i++) {
			if (final_values[i][1] > 1e-6) {
				fprintf(output, pk_format, final_values[i][0], final_values[i][1], final_values[i][2], final_values[i][3], final_kBins[i], final_values[i][4], final_values[i][5]);
			} else {
				printf("monopole 0 in fold() ? \n");
				exit(0);
			}
		}	
	}	
	
	fclose(unfolded_file);
	fclose(folded_file);
	fclose(output);
	free2D_double(&folded_values, len_folded);
	free2D_double(&unfolded_values, len_unfolded);
	return toReturn;
}



//THIS PROBABLY DOESNT WORK
/* 
int combine_multiple(char paths_unfolded[], char paths_folded[], char outFiles[], double kmin, double kmax, double k_fold, int simNo, int first_index) {

	char pu[100];
	char pf[100];
	char out[100];
	FoldingInformation allInfo[simNo];
	FoldingInformation thisInfo;
			
	//collect info for all sims - where the transition from unfolded to folded happens
	for (int i = first_index; i < simNo + first_index; i++) {

		sprintf(pu, paths_unfolded, i);
		sprintf(pf, paths_folded, i);
		sprintf(out, outFiles, i);		
		//get indices, dont print
		//printf("%s\n", pu);
		thisInfo = combine_folded(pu, pf, out, kmin, kmax, k_fold, false);
		allInfo[i - first_index] = thisInfo;	
	}
	
	//sum how many sims jump from [ind] index in the unfolded file
	//[bin][0 - how many counts of this starting index, 1 - which index it connects to, 2 - kmin ind, 3 - foldover]
	int index_popularity[100][4];
	int ind;

	for (int i = 0; i < simNo; i++) {
		thisInfo = allInfo[i];
		ind = thisInfo.unfold_last_ind;
		index_popularity[ind][0] += 1;
		index_popularity[ind][1] = thisInfo.fold_first_ind;
		index_popularity[ind][2] = thisInfo.min_k_ind;		
		index_popularity[ind][3] = thisInfo.foldover_index;	
	}
	
	//find which [ind] index has greatest number of unfolded files ending there
	int ind_end_unfold, ind_start_fold, maxCount, maxAt, min_k_ind, foldover_index;
	maxCount = 0;
	for (int i = 0; i < 100; i++) {
		if (index_popularity[i][0] > maxCount) {
			maxCount = index_popularity[i][0];
			maxAt = i;
		}	
	}

	
	//set final variables
	ind_end_unfold = maxAt;
	ind_start_fold = index_popularity[maxAt][1];
	min_k_ind = index_popularity[maxAt][2];
	foldover_index = index_popularity[maxAt][3];
	
	FILE *final_output;
	//print the necessary files	
	for (int sim = 0; sim < simNo; sim++) {	
		thisInfo = allInfo[sim];
		sprintf(out, outFiles, sim + first_index);		
		final_output = fopen(out, "w");	
		double final_values[200][4];
		for (int i = min_k_ind; i < ind_end_unfold; i++) {
			for (int k = 0; k < 4; k++) final_values[i - min_k_ind][k] = thisInfo.unfolded_values[i][k];
		}	
		
		//printf("ind start fold: %d, foldover: %d, ind_end_unfold: %d\n", ind_start_fold, foldover_index, ind_end_unfold);		

		for (int i = ind_start_fold; i < foldover_index; i++) {
			for (int k = 0; k < 4; k++) {			
				final_values[ind_end_unfold - min_k_ind + i - ind_start_fold][k] = thisInfo.folded_values[i][k];
			}	
		}	
		
		for (int i = 0; i < 200; i++) {
			if (final_values[i][1] > 0.000001) {
				fprintf(final_output, "%.5lf \t %.5lf \t %.5lf \t %d\n", final_values[i][0], final_values[i][1], final_values[i][2], (int) final_values[i][3]);
			}
		}	
		fclose(final_output);	
	}
	return 0;
}

*/
