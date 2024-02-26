double muOrderZero(double ks){
    // limits were established with signa = 2.*3./sqrt(2.), limits should scale propto velDispersion/(2.*3./sqrt(2.))
    if(ks > 0.0914289)  return 0.5*pow(pi, 0.5)*gsl_sf_erf(ks)/ks;

    else{
        return 1.0 - 0.333333*pow(ks, 2.) + 0.1*pow(ks, 4.) - 0.0238095*pow(ks, 6.);
    }
}


double muOrderTwo(double ks){
    if(ks > 0.15727469)  return pow(4.*ks*ks*ks, -1.)*(pow(pi, 0.5)*gsl_sf_erf(ks) - 2.*ks*exp(-1.*ks*ks));
    
    else{
        return 1./3. - ks*ks/5. + pow(ks, 4.)/14.;
    }
}


double muOrderFour(double ks){
    if(ks>0.2418305)  return pow(8.*pow(ks, 5.), -1.)*(3.*pow(pi, 0.5)*gsl_sf_erf(ks) -6.*ks*exp(-1.*ks*ks) - 4.*pow(ks, 3.)*exp(-ks*ks));

    else{
        return 1./5. - ks*ks/7. + pow(ks, 4.)/18.;
    }
}


double muOrderSix(double ks){
    if(ks>0.335168)  return pow(16.*pow(ks, 7.), -1.)*(15.*pow(pi, 0.5)*gsl_sf_erf(ks) - 2.*ks*exp(-ks*ks)*(15. + 10.*ks*ks + 4.*pow(ks, 4.)));

    else{
         return 1./7. - ks*ks/9. + pow(ks, 4.)/22.;
    }
}


double muOrderEight(double ks){
    if(ks>0.4326645)  return (105.*pow(pi, 0.5)*gsl_sf_erf(ks)*pow(32.*pow(ks, 9.), -1.) - 2.*ks*exp(-ks*ks)*pow(32.*pow(ks, 9.), -1.)*(8.*pow(ks, 6.) + 28.*pow(ks, 4.) + 70.*pow(ks, 2.) + 105.));

    else{ 
        // Numerical result diverges for ks < 0.2, replace by Taylor expansion at ksigma=0.5
         return 1./9. - ks*ks/11. + pow(ks, 4.)/26.;
    }
}



// Multipoles for the Kaiser-Gauss redshift space distortion model.

double kaiserGauss_Monofactor(double ks, double beta){
    return muOrderZero(ks) + 2.*beta*muOrderTwo(ks) + beta*beta*muOrderFour(ks);
}


double kaiserGauss_Quadfactor(double ks, double beta){
    return (5./2.)*(-1.*muOrderZero(ks) + muOrderTwo(ks)*(3. - 2.*beta) + muOrderFour(ks)*(-beta*beta + 6.*beta) + 3.*beta*beta*muOrderSix(ks));
}

double Kaiser_Monofactor(double beta) {
	return 1.0 + (2.0/3.0)*beta + (1.0/5.0)*pow(beta, 2.0);
}

double Kaiser_Quadfactor(double beta) {
	return (4.0/3.0)*beta + (4.0/7.0)*pow(beta, 2.0);
}

double Kaiser_Hexfactor(double beta) {
	return (8.0/35.0)*pow(beta, 2.0);
}

