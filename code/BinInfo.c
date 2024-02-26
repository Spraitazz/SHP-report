

//a single struct to pass around binning stuff
BinInfo *prep_bins(double xmin, double xmax, int bins, int bin_type) {
	BinInfo* toReturn = malloc(sizeof(*toReturn));
	if (toReturn == NULL) {
		printf("error in prep_bins()\n");
		exit(0);
	} else {
		toReturn->bins = bins;
		if (bin_type == LOG_BIN) {
			double logxmin = 0.0;
			if (xmin < 1e-10) {
				logxmin = xmin;
			} else {
				logxmin = log(xmin);
			}
			double logxmax = log(xmax);
			double logdx = (logxmax - logxmin) / (double)(bins-1);
			toReturn->xmin = logxmin;
			toReturn->xmax = logxmax;
			toReturn->dx = logdx;
		} else if (bin_type == NORMAL_BIN) {
			double dx = (xmax - xmin) / (double)(bins-1);
			toReturn->xmin = xmin;
			toReturn->xmax = xmax;
			toReturn->dx = dx;
		} else {
			printf("bin types can only be NORMAL_BIN and LOG_BIN in prep_bins() \n");
			exit(0);
		}
		toReturn->bin_type = bin_type;
	}
	return toReturn;
}

/*

//ACTUALLY NEED: more general object Bins that is just an array of x values
// and always use that, BinInfo can be used to prepare one Bins

BinInfo *copy_bins(Spline *s) {

}
*/
int x_to_bin(BinInfo *binInfo, double x) {
	int index = 0;
	if (binInfo->bin_type == LOG_BIN) {
		index = (int) floor((log(x) - binInfo->xmin) / binInfo->dx);
	} else {
		index = (int) floor((x - binInfo->xmin) / binInfo->dx);
	}
	if (index < 0 || index >= binInfo->bins) {
		printf("x_to_bin index %d, bins: %d \n", index, binInfo->bins);
		exit(0);
	}
	return index;
}

double bin_to_x(BinInfo *binInfo, int index) {
	double value = (double)index * binInfo->dx + binInfo->xmin;
	if (binInfo->bin_type == LOG_BIN) value = exp(value);	
	return value;
}

int print_binInfo(BinInfo *binInfo) {
	if (binInfo->bin_type == NORMAL_BIN) {
		printf("%d evenly spaced bins from %le to %le\n", binInfo->bins, binInfo->xmin, binInfo->xmax); 
	} else {
		printf("%d evenly LOG-spaced bins from %le to %le\n", binInfo->bins, exp(binInfo->xmin), exp(binInfo->xmax)); 
	}
	return 0;
}

