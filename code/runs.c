int haloes_measure_Pk(Catalogue *catalogue, char outFile[], double shotnoise, double foldfactor, int grid_func, double multip_fact) {
	
	bool folded = false;
	Catalogue *particles_saved;
	
	if (foldfactor < 1.0) {
		printf("what are you doing? foldfactor: %lf\n", foldfactor);
		exit(0);
	} else if (foldfactor > 1.0) {
		Jenkins_foldfactor = foldfactor;
		folded = true;
		//particles_saved = Jenkins_fold_volume(catalogue);		
	}
	
	printf("Measuring P(k) of %d particles with %lf factor for overdensity. Subtracting shotnoise equal to %lf\n", catalogue->particle_no, multip_fact, shotnoise);
	
	
	//allocate grid memory
	initGrid();	
	
	init_spectrum_storage(spectrum_size);

	populateGrid(catalogue, grid_func);	

	overdensity_allocateMemory();	
	gridIntoOverdensity(catalogue->particle_no);	
	overdensity_pk(multip_fact);

	pk_to_file(outFile, shotnoise);	

	//if (folded)	Jenkins_restore(catalogue, particles_saved);

	free_spectrum_storage();
	freeGrid();
	overdensity_freeMemory();
	return 0;
}

int catalogues_measure_Pk(char inFiles[], char outFiles[], char format[], int cat_no, int startInd, double foldfactor, int grid_func, double multip_fact, Parameters *params) {
	
	initGrid();	
	overdensity_allocateMemory();	
	init_spectrum_storage(spectrum_size);
	
	char cat_in[200];
	char Pk_out[200];
	Catalogue *catalogue = NULL;
	Catalogue *particles_saved = NULL;	
	
	for (int i = 0; i < cat_no; i++) {
		
		sprintf(cat_in, inFiles, i+startInd);
		catalogue = input_catalogue_file(cat_in, 0, format);
		toRedshift(catalogue, params);
	
		bool folded = false;	
		if (foldfactor < 1.0) {
			printf("what are you doing? foldfactor: %lf\n", foldfactor);
			exit(0);
		} else if (foldfactor > 1.0) {
			Jenkins_foldfactor = foldfactor;
			folded = true;
			//particles_saved = Jenkins_fold_volume(catalogue);		
		}
	
		printf("Measuring P(k) of %d particles with %lf factor for overdensity.\n", catalogue->particle_no, multip_fact);
			

		populateGrid(catalogue, grid_func);	

		gridIntoOverdensity(catalogue->particle_no);	
		overdensity_pk(multip_fact);

		sprintf(Pk_out, outFiles, i+startInd);	
		pk_to_file(Pk_out, 0.0);	

		//if (folded)	Jenkins_restore(catalogue, particles_saved);
		
		clearPk();
		overdensity_clearMemory();
		clearGrid();
		
		free_catalogue(catalogue);		
	
	}
	
	free_spectrum_storage();
	freeGrid();
	overdensity_freeMemory();	
	return 0;
}

/*
Foldings* measure_Pk(Foldings* foldings, double kmin, double kmax, char outFile[]) {

	//printf("Pk folding from kmin: %lf to kmax: %lf \n", kmin, kmax);
	
	char tmp_unfolded[100], tmp_folded[100], tmp_combined[100];
	FoldingInformation* foldingInfo;
	
	//CHECK FOR KMIN - IF TOO BIG, START STRAIGHT FROM FOLDED
	
	sprintf(tmp_unfolded, "%s/data/tmp/tmp_uf.dat", home_directory);
	sprintf(tmp_folded, "%s/data/tmp/tmp_f.dat", home_directory);
	sprintf(tmp_combined, "%s/data/tmp/tmp_comb.dat", home_directory);
	
	double foldFactor = 2.0;
	Foldings *toReturn = malloc(sizeof(*toReturn));
	bool max_reached = false;
	
	if (foldings == NULL) {
	
		int folds = 0;
		int tmp_foldPoints[50][4];
		
		//haloes_measure_Pk(tmp_unfolded, 1.0, NGP);
		//haloes_measure_Pk(tmp_folded, 2.0, NGP);
		foldingInfo = combine_folded(tmp_unfolded, tmp_folded, pk_format, tmp_combined, kmin, kmax, -1.0, NULL, true);
		max_reached = foldingInfo->max_reached;
		
		tmp_foldPoints[folds][0] = foldingInfo->uf_first_ind;
		tmp_foldPoints[folds][1] = foldingInfo->uf_last_ind;
		tmp_foldPoints[folds][2] = foldingInfo->f_first_ind;
		tmp_foldPoints[folds][3] = foldingInfo->f_last_ind;
		
		folds += 1;
		remove(tmp_unfolded);
		remove(tmp_folded);
		rename(tmp_combined, tmp_unfolded);	
		free(foldingInfo);		
		
		while (!max_reached) {			
			
			foldFactor += 2.0;		
			//printf("foldfactor: %lf \n", foldFactor);
			
			//haloes_measure_Pk(tmp_folded, foldFactor, NGP);
			foldingInfo = combine_folded(tmp_unfolded, tmp_folded, pk_format, tmp_combined, kmin, kmax, -1.0, NULL, true);
			max_reached = foldingInfo->max_reached;
			
			tmp_foldPoints[folds][0] = foldingInfo->uf_first_ind;
			tmp_foldPoints[folds][1] = foldingInfo->uf_last_ind;
			tmp_foldPoints[folds][2] = foldingInfo->f_first_ind;
			tmp_foldPoints[folds][3] = foldingInfo->f_last_ind;
			
			folds += 1;
			remove(tmp_unfolded);
			remove(tmp_folded);
			rename(tmp_combined, tmp_unfolded);
			free(foldingInfo);	
		
		}		
		
		rename(tmp_unfolded, outFile);
		
		int** foldPoints;
		calloc2D_int(&foldPoints, folds, 4);
		for (int i = 0; i < folds; i++) {
			//printf("connection %d. uf first: %d, uf last: %d, f first: %d, f last: %d \n", i, tmp_foldPoints[i][0], tmp_foldPoints[i][1], tmp_foldPoints[i][2], tmp_foldPoints[i][3]);
			foldPoints[i][0] = tmp_foldPoints[i][0];
			foldPoints[i][1] = tmp_foldPoints[i][1];
			foldPoints[i][2] = tmp_foldPoints[i][2];
			foldPoints[i][3] = tmp_foldPoints[i][3];
		}			
		
		toReturn->foldPoints = foldPoints;
		toReturn->folds = folds;		
	
	} else {	
		
		//haloes_measure_Pk(tmp_unfolded, 1.0, NGP);
		bool success = true;
	
		for (int i = 0; i < foldings->folds; i++) {
		
			
			//int foldInfo[4];			
			//foldInfo[0] = foldings->foldPoints[i][0]; 
			//foldInfo[1] = foldings->foldPoints[i][1]; 
			//foldInfo[2] = foldings->foldPoints[i][2]; 
			//foldInfo[3] = foldings->foldPoints[i][3]; 
			
			int* foldInfo = foldings->foldPoints[i];			
			
			//haloes_measure_Pk(tmp_folded, foldFactor, NGP);
			//printf("combining, foldfactor: %lf, uf end: %d, f start: %d \n", foldFactor, foldInfo[0], foldInfo[1]);
			foldingInfo = combine_folded(tmp_unfolded, tmp_folded, pk_format, tmp_combined, kmin, kmax, -1.0, foldInfo, true);
			
			if(!foldingInfo->success) {
				remove(tmp_unfolded);
				remove(tmp_folded);
				remove(tmp_combined);
				return foldings;
			}
			
			free(foldingInfo);
			//printf("combined, foldfactor: %lf \n", foldFactor);
			
			remove(tmp_unfolded);	
			remove(tmp_folded);
			rename(tmp_combined, tmp_unfolded);
			foldFactor += 2.0;	
		
		}
		
		toReturn = foldings;		
		rename(tmp_unfolded, outFile);
	
	}
	return toReturn;
}
*/





