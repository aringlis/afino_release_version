
"""
Main analysis module for Bayesian MCMC Fourier spectrum fitting
"""

import numpy as np
import copy
import os
import sys
from matplotlib import pyplot as plt
from afino_series import AfinoSeries
import rnspectralmodels3
import afino_model_fitting
import pickle

def main_analysis(ts,model='single_power_law_with_constant',low_frequency_cutoff=None,overwrite_gauss_bounds = None, overwrite_extra_gauss_bounds = None):



    #get the Fourier spectrum out of the timeseries object
    nt = ts.SampleTimes.nt
    dt = ts.SampleTimes.dt
    iobs = ts.PowerSpectrum.Npower #change this to .npower for normalized version. Use this instead, could be an issue with priors.
    this = ([ts.PowerSpectrum.frequencies.positive, iobs],)
    frequencies = ts.PowerSpectrum.frequencies.positive

    if low_frequency_cutoff:
        mask = [frequencies < low_frequency_cutoff]
        frequencies = frequencies[mask]
        iobs = iobs[mask]
    
    #---------------------------------------
    #use fitting methods from scipy optimize
    best_lnlike = -999999.0

    # --------------------------------------------------------
    if model == 'single_power_law_with_constant':

        #don't have a good idea of starting guess parameters so randomize these and do multiple trials to cover more parameter space
        for i in range(0,20):

            #get randomized initial guess params to input into fit
            guess = randomize_initial_guess()

            #try 3 different fitting algorithms to ensure we maximize the likelihood
            for method in ['L-BFGS-B','TNC','SLSQP']:
                res = afino_model_fitting.go_plaw(frequencies,iobs,rnspectralmodels3.power_law_with_constant,guess[0:3],method)
                param_vals = res['x']
                jack_lnlike = afino_model_fitting.lnlike(param_vals,frequencies,iobs,rnspectralmodels3.power_law_with_constant)
            

                if jack_lnlike > best_lnlike:
                    best_lnlike = copy.deepcopy(jack_lnlike)
                    best_param_vals = copy.deepcopy(param_vals)

        #calculate BIC and store best fit power spectrum and params
        jack_bic = afino_model_fitting.BIC(2,best_param_vals,frequencies,iobs,rnspectralmodels3.power_law_with_constant,len(iobs))
        best_fit_power_spectrum = rnspectralmodels3.power_law_with_constant(best_param_vals,frequencies)
        best_fit_params = best_param_vals
            

        
    #-------------------------------------    
    if model == 'splwc_AddNormalBump2':
        
        for i in range(0,20):

            #get randomized initial guess params to input into fit
            guess = randomize_initial_guess()

            #try 3 different fitting algorithms to ensure we maximize the likelihood
            for method in ['L-BFGS-B','TNC','SLSQP']:
                res = afino_model_fitting.go_gauss(frequencies,iobs,rnspectralmodels3.splwc_AddNormalBump2,guess,method,overwrite_gauss_bounds = overwrite_gauss_bounds)
                param_vals = res['x']
                jack_lnlike = afino_model_fitting.lnlike(param_vals,frequencies,iobs,rnspectralmodels3.splwc_AddNormalBump2)

                if jack_lnlike > best_lnlike:
                    best_lnlike = copy.deepcopy(jack_lnlike)
                    best_param_vals = copy.deepcopy(param_vals)

        #calculate BIC and store best fit power spectrum and params
        jack_bic = afino_model_fitting.BIC(5,best_param_vals,frequencies,iobs,rnspectralmodels3.splwc_AddNormalBump2,len(iobs))
        best_fit_power_spectrum = rnspectralmodels3.splwc_AddNormalBump2(best_param_vals,frequencies)
        best_fit_params = best_param_vals
          
        

    # --------------------------------------------------------
    if model == 'broken_power_law_with_constant':

        #don't have a good idea of starting guess parameters so randomize these and do multiple trials to cover more parameter space
        for i in range(0,20):

            #get randomized initial guess params to input into fit
            guess = randomize_initial_guess_bpow()

            #try 3 different fitting algorithms to ensure we maximize the likelihood
            for method in ['L-BFGS-B','TNC','SLSQP']:
                res = afino_model_fitting.go_bpow(frequencies,iobs,rnspectralmodels3.broken_power_law_with_constant,guess,method)
                param_vals = res['x']
                jack_lnlike = afino_model_fitting.lnlike(param_vals,frequencies,iobs,rnspectralmodels3.broken_power_law_with_constant)

                if jack_lnlike > best_lnlike:
                    best_lnlike = copy.deepcopy(jack_lnlike)
                    best_param_vals = copy.deepcopy(param_vals)

        #calculate BIC and store best fit power spectrum and params
        jack_bic = afino_model_fitting.BIC(4,best_param_vals,frequencies,iobs,rnspectralmodels3.broken_power_law_with_constant,len(iobs))
        best_fit_power_spectrum = rnspectralmodels3.broken_power_law_with_constant(best_param_vals,frequencies)
        best_fit_params = best_param_vals

                    

    # add additional models as needed here
    #-------------------------------------
    
    #------------------------------------------------------------
    if model == 'splwc_AddNormalBump2_plus_extra_bump':

        #don't have a good idea of starting guess parameters so randomize these and do multiple trials to cover more parameter space
        for i in range(0,20):

            #get randomized initial guess params to input into fit
            guess = randomize_initial_guess_plus_extra_bump()

            #try 3 different fitting algorithms to ensure we maximize the likelihood
            for method in ['L-BFGS-B','TNC','SLSQP']:
                res = afino_model_fitting.go_gauss_plus_extra_bump(frequencies,iobs,rnspectralmodels3.splwc_AddNormalBump2_plus_extra_bump,guess,method, overwrite_extra_gauss_bounds = overwrite_extra_gauss_bounds)
                param_vals = res['x']
                jack_lnlike = afino_model_fitting.lnlike(param_vals,frequencies,iobs,rnspectralmodels3.splwc_AddNormalBump2_plus_extra_bump)

                if jack_lnlike > best_lnlike:
                    best_lnlike = copy.deepcopy(jack_lnlike)
                    best_param_vals = copy.deepcopy(param_vals)

        #calculate BIC and store best fit power spectrum and params
        jack_bic = afino_model_fitting.BIC(8,best_param_vals,frequencies,iobs,rnspectralmodels3.splwc_AddNormalBump2_plus_extra_bump,len(iobs))
        best_fit_power_spectrum = rnspectralmodels3.splwc_AddNormalBump2_plus_extra_bump(best_param_vals,frequencies)
        best_fit_params = best_param_vals



    #----------------------------------------
            
    #calculate the uncertainty on the parameters using the inverse of the Hessian
   # pickle.dump([frequencies,iobs,model],open('file_for_hessian.pkl','wb'))
   # h = numdifftools.Hessian(lnlike_for_hessian)
   # hess = h(best_fit_params)
   # hess_inv = np.linalg.inv(hess)
   # uncertainties = np.sqrt(np.abs(np.diag(hess_inv)))


    #also want to find the goodness of fit for the best model
    smr = afino_model_fitting.rhoj(iobs,best_fit_power_spectrum)
    if model == 'single_power_law_with_constant':
        deg_free = len(iobs) - 2
        rchi2 = afino_model_fitting.rchi2(1, deg_free, smr)
        prob = afino_model_fitting.prob_this_rchi2_or_larger(rchi2, 1, deg_free)
    elif model == 'broken_power_law_with_constant':
        deg_free = len(iobs) - 4
        rchi2 = afino_model_fitting.rchi2(1, deg_free, smr)
        prob = afino_model_fitting.prob_this_rchi2_or_larger(rchi2, 1, deg_free)
    elif model == 'splwc_AddNormalBump2_plus_extra_bump':
        deg_free = len(iobs) - 8
        rchi2 = afino_model_fitting.rchi2(1, deg_free, smr)
        prob = afino_model_fitting.prob_this_rchi2_or_larger(rchi2, 1, deg_free)
    else:
        deg_free = len(iobs) - 5
        rchi2 = afino_model_fitting.rchi2(1,deg_free,smr)
        prob = afino_model_fitting.prob_this_rchi2_or_larger(rchi2, 1, deg_free)


    #now have the best fit, likelihood and BIC value. Return all these in a dictionary
    #will save combined fit results (both models) to a pickle file in the next level up.
    fitresults = {}
    fitresults['lnlike'] = best_lnlike #likelihoods[selection_index]
    fitresults['model'] = model
    fitresults['BIC'] = jack_bic
    fitresults['best_fit_power_spectrum'] = best_fit_power_spectrum
    fitresults['frequencies'] = frequencies
    fitresults['power'] = iobs
    fitresults['params'] = best_fit_params
    fitresults['rchi2'] = rchi2
    fitresults['probability'] = prob
    
    return fitresults
            
    
        
 

def randomize_initial_guess():
    
    #randomize width between 0.05 and 0.25
    gauss_width = (np.random.random(1) / 5.) + 0.05
    gauss_position = np.random.random(1) * (-4) - 1.5
    gauss_amp = np.random.random(1) * (-14)
    plaw_amp = (np.random.random(1) *20) - 10
    plaw_index = np.random.random(1) * (-5) - 1
    background = np.random.random(1) * (-20)

    initial_guess = [plaw_amp[0], plaw_index[0], background[0], gauss_amp[0], gauss_position[0], gauss_width[0]]

    return initial_guess
    

def randomize_initial_guess_plus_extra_bump():
    
    #randomize width between 0.05 and 0.25
    gauss_width = (np.random.random(1) / 5.) + 0.05
    gauss_position = np.random.random(1) * (-4) - 1.5
    gauss_amp = np.random.random(1) * (-14)
    plaw_amp = (np.random.random(1) *20) - 10
    plaw_index = np.random.random(1) * (-5) - 1
    background = np.random.random(1) * (-20)
    
    extra_gauss_amp = np.random.random(1) * (-14)
    extra_gauss_position = (np.random.random(1) / 5.) -3.1
    extra_gauss_width = (np.random.random(1) / 10.0) + 0.05
    
    initial_guess = [plaw_amp[0], plaw_index[0], background[0], gauss_amp[0], gauss_position[0], gauss_width[0], extra_gauss_amp[0], extra_gauss_position[0], extra_gauss_width[0]]

    return initial_guess


def randomize_initial_guess_bpow():
     
    plaw_amp = (np.random.random(1) *20) - 10
    plaw_index = np.random.random(1) * (-5) - 1
   # plaw_index2 = np.random.random(1) * (-5) - 1
    break_position = np.random.random(1) * (-4) - 1.5
    background = np.random.random(1) * (-20)

    initial_guess = [plaw_amp[0], plaw_index[0], np.exp(break_position[0]), plaw_index[0], background[0]]

    return initial_guess


