


//Kaiser-Gauss
/*
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
*/

double muOrderZero(double ks){
    if(ks > 0.04)  return pow(2., 0.5)*atan(ks*pow(2., -0.5))/ks;

    else{
      return 1.0 - 0.166667*pow(ks, 2.) + 0.05*pow(ks, 4.) - 0.0178571*pow(ks, 6.);
    }
}

double muOrderTwo(double ks){
      if(ks > 0.051675)  return 2./pow(ks,2.) - 2.*sqrt(2.)*atan(ks*pow(2., -0.5))/pow(ks, 3.); 

      else{
        return 1./3. - 0.1*pow(ks, 2.) + 0.0357143*pow(ks, 4.) - 0.0138889*pow(ks, 6.);
      }
}

double muOrderFour(double ks){
      if(ks > 0.0512935)  return -4./pow(ks, 4.) +(2./3.)/pow(ks, 2.) + 4.*sqrt(2.)*atan(ks/sqrt(2.))/pow(ks, 5.); 

      else{
        // First term in taylor expansion of M_n is 1/(n+1) i.e. Mn(sigma=0) = 1/(n+1)
        return 0.2 - 0.0714286*pow(ks, 2.) + 0.0277778*pow(ks, 4.) - 0.0113636*pow(ks, 6.);
      }
}

double muOrderSix(double ks){
      if(ks > 0.039329)  return 8./pow(ks,6.) - (4./3.)*pow(ks, -4.) + 0.4*pow(ks, -2.) -8.*sqrt(2.)*atan(ks/sqrt(2.))/pow(ks, 7.);

      else{
         return 0.142857 - 0.0555556*pow(ks, 2.) + 0.0227273*pow(ks, 4.) - 0.00961538*pow(ks, 6.);
      }
}

double muOrderEight(double ks){
    if(ks>0.181839)  return -16.*pow(ks,-8.) + (8./3.)*pow(ks, -6.) - 0.8*pow(ks, -4.) + (2./7.)*pow(ks, -2.) + 16.*sqrt(2.)*atan(ks/sqrt(2.))/pow(ks, 9.);

    else{
      return 1./9. - 0.0454545*ks*ks + 0.0192308*pow(ks, 4.) - 0.00833333*pow(ks, 6.);
    }
}

double muOrderTen(double ks){
    if(ks>0.14)  return (2./315.)*(35./pow(ks,2.) -90./pow(ks,4.) +252./pow(ks,6.) -840./pow(ks,8.) +5040./pow(ks,10.) -5040.*sqrt(2.)*atan(ks/sqrt(2.))/pow(ks,11.));

    else{
        return 1./11. - ks*ks/26. + pow(ks, 4.)/60. - pow(ks, 6.)/136. + pow(ks, 8.)/304.;
    }
}

// Multipoles for the Kaiser-Lorentz redshift space distortion model.
double KaiserLorentz_Monofactor(double ks, double beta){
    return muOrderZero(ks) + 2.*beta*muOrderTwo(ks) + beta*beta*muOrderFour(ks);
}

double KaiserLorentz_Quadfactor(double ks, double beta){
    return (5./2.)*(-1.*muOrderZero(ks) + muOrderTwo(ks)*(3. - 2.*beta) + muOrderFour(ks)*(-beta*beta + 6.*beta) + 3.*beta*beta*muOrderSix(ks));
}

double KaiserLorentz_Hexfactor(double ks, double beta){
    return (9./8.)*(35.*beta*beta*muOrderEight(ks)  + 10.*beta*(7. -3.*beta)*muOrderSix(ks) + (35. - 60.*beta + 3.*beta*beta)*muOrderFour(ks) + 6.*(beta - 5.)*muOrderTwo(ks) + 3.*muOrderZero(ks));
}

//Mike
double kaiserGauss_Monofactor(double ks, double beta){
    return muOrderZero(ks) + 2.*beta*muOrderTwo(ks) + beta*beta*muOrderFour(ks);
}

//Mike
double kaiserGauss_Quadfactor(double ks, double beta){
    return (5./2.)*(-1.*muOrderZero(ks) + muOrderTwo(ks)*(3. - 2.*beta) + muOrderFour(ks)*(-beta*beta + 6.*beta) + 3.*beta*beta*muOrderSix(ks));
}

double kaiserGauss_Hexfactor(double ks, double beta){
    return (9./8.)*(35.*beta*beta*muOrderEight(ks)  -30.*beta*beta*muOrderSix(ks) + 3.*beta*beta*muOrderFour(ks) + 70.*beta*muOrderSix(ks) - 60.*beta*muOrderFour(ks)  + 6.*beta*muOrderTwo(ks) + 35.*muOrderFour(ks) - 30.*muOrderTwo(ks) + 3.*muOrderZero(ks));
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




