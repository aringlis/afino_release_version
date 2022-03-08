
"""
This routine initiates the AFINO analysis method for a particular chosen model.
The fitting is done via SciPy, and the results, including best fit and associated
BIC value, are stored in a dictionary and returned.
"""

import numpy as np
import copy
from afino import afino_spectral_models
from afino import afino_model_fitting

def main_analysis(ts, model='pow_const', low_frequency_cutoff=None,
                      overwrite_gauss_bounds = None, overwrite_extra_gauss_bounds = None):



    # get the Fourier spectrum out of the timeseries object
    nt = ts.SampleTimes.nt
    dt = ts.SampleTimes.dt
    iobs = ts.PowerSpectrum.Npower 
    this = ([ts.PowerSpectrum.frequencies.positive, iobs],)
    frequencies = ts.PowerSpectrum.frequencies.positive

    # if a low frequency cutoff is set, only analyse the spectrum at frequencies above it
    if low_frequency_cutoff:
        mask = [frequencies < low_frequency_cutoff]
        frequencies = frequencies[tuple(mask)]
        iobs = iobs[tuple(mask)]
    
    #---------------------------------------
    #use fitting methods from scipy optimize
    best_lnlike = -999999.0

    # --------------------------------------------------------


    #don't have a good idea of starting guess parameters so randomize these and do multiple trials to cover more parameter space
    for i in range(0,20):

        #get randomized initial guess params to input into fit
        guess = randomize_initial_guess(model = model)

        #try 3 different fitting algorithms to ensure we maximize the likelihood
        for method in ['L-BFGS-B','SLSQP']:
            if model == 'pow_const':
                res = afino_model_fitting.go_plaw(frequencies,iobs,afino_spectral_models.pow_const,guess,method)
                param_vals = res['x']
                jack_lnlike = afino_model_fitting.lnlike(param_vals,frequencies,iobs,afino_spectral_models.pow_const)
                if jack_lnlike > best_lnlike:
                    best_lnlike = copy.deepcopy(jack_lnlike)
                    best_param_vals = copy.deepcopy(param_vals)               

                
            elif model == 'pow_const_gauss':
                res = afino_model_fitting.go_gauss(frequencies,iobs,afino_spectral_models.pow_const_gauss,guess,method,overwrite_gauss_bounds = overwrite_gauss_bounds)
                param_vals = res['x']
                jack_lnlike = afino_model_fitting.lnlike(param_vals,frequencies,iobs,afino_spectral_models.pow_const_gauss)
                if jack_lnlike > best_lnlike:
                    best_lnlike = copy.deepcopy(jack_lnlike)
                    best_param_vals = copy.deepcopy(param_vals)


            elif model == 'bpow_const':
                res = afino_model_fitting.go_bpow(frequencies,iobs,afino_spectral_models.bpow_const,guess,method)
                param_vals = res['x']
                jack_lnlike = afino_model_fitting.lnlike(param_vals,frequencies,iobs,afino_spectral_models.bpow_const)
                if jack_lnlike > best_lnlike:
                    best_lnlike = copy.deepcopy(jack_lnlike)
                    best_param_vals = copy.deepcopy(param_vals)


            elif model == 'pow_const_2gauss':
                res = afino_model_fitting.go_gauss_plus_extra_bump(frequencies,iobs,afino_spectral_models.pow_const_2gauss,guess,method,
                                                                   overwrite_extra_gauss_bounds = overwrite_extra_gauss_bounds)
                param_vals = res['x']
                jack_lnlike = afino_model_fitting.lnlike(param_vals,frequencies,iobs,afino_spectral_models.pow_const_2gauss)
                if jack_lnlike > best_lnlike:
                    best_lnlike = copy.deepcopy(jack_lnlike)
                    best_param_vals = copy.deepcopy(param_vals)
                
            else:
                raise ValueError
            

    #calculate BIC and store best fit power spectrum and params
    #also want to find the goodness of fit for the best model
    if model == 'pow_const':
        jack_bic = afino_model_fitting.BIC(2,best_param_vals,frequencies,iobs,afino_spectral_models.pow_const,len(iobs))
        best_fit_power_spectrum = afino_spectral_models.pow_const(best_param_vals,frequencies)
        smr = afino_model_fitting.rhoj(iobs,best_fit_power_spectrum)
        deg_free = len(iobs) - 2
        rchi2 = afino_model_fitting.rchi2(1, deg_free, smr)
        prob = afino_model_fitting.prob_this_rchi2_or_larger(rchi2, 1, deg_free)
    elif model == 'pow_const_gauss':
        jack_bic = afino_model_fitting.BIC(5,best_param_vals,frequencies,iobs,afino_spectral_models.pow_const_gauss,len(iobs))
        best_fit_power_spectrum = afino_spectral_models.pow_const_gauss(best_param_vals,frequencies)
        smr = afino_model_fitting.rhoj(iobs,best_fit_power_spectrum)
        deg_free = len(iobs) - 5
        rchi2 = afino_model_fitting.rchi2(1,deg_free,smr)
        prob = afino_model_fitting.prob_this_rchi2_or_larger(rchi2, 1, deg_free)
    elif model == 'bpow_const':
        jack_bic = afino_model_fitting.BIC(4,best_param_vals,frequencies,iobs,afino_spectral_models.bpow_const,len(iobs))
        best_fit_power_spectrum = afino_spectral_models.bpow_const(best_param_vals,frequencies)
        smr = afino_model_fitting.rhoj(iobs,best_fit_power_spectrum)
        deg_free = len(iobs) - 4
        rchi2 = afino_model_fitting.rchi2(1, deg_free, smr)
        prob = afino_model_fitting.prob_this_rchi2_or_larger(rchi2, 1, deg_free)
    elif model == 'pow_const_2gauss':
        jack_bic = afino_model_fitting.BIC(8,best_param_vals,frequencies,iobs,afino_spectral_models.pow_const_2gauss,len(iobs))
        best_fit_power_spectrum = afino_spectral_models.pow_const_2gauss(best_param_vals,frequencies)
        smr = afino_model_fitting.rhoj(iobs,best_fit_power_spectrum)
        deg_free = len(iobs) - 8
        rchi2 = afino_model_fitting.rchi2(1, deg_free, smr)
        prob = afino_model_fitting.prob_this_rchi2_or_larger(rchi2, 1, deg_free)
    best_fit_params = best_param_vals


    #----------------------------------------
    

    #now have the best fit, likelihood and BIC value. Return all these in a dictionary
    #will save combined fit results (both models) to a pickle file in the next level up.
    fitresults = {}
    fitresults['lnlike'] = best_lnlike 
    fitresults['model'] = model
    fitresults['BIC'] = jack_bic
    fitresults['best_fit_power_spectrum'] = best_fit_power_spectrum
    fitresults['frequencies'] = frequencies
    fitresults['power'] = iobs
    fitresults['params'] = best_fit_params
    fitresults['rchi2'] = rchi2
    fitresults['probability'] = prob
    
    return fitresults
            
    
        
 

def randomize_initial_guess(model = 'pow_const'):
    
    # randomize initial guesses for all model starting parameters
    gauss_width = (np.random.random(1) / 5.) + 0.05
    gauss_position = np.random.random(1) * (-4) - 1.5
    gauss_amp = np.random.random(1) * (-14)
    plaw_amp = (np.random.random(1) *20) - 10
    plaw_index = np.random.random(1) * (-5) - 1
    background = np.random.random(1) * (-20)
    
    second_gauss_amp = np.random.random(1) * (-14)
    second_gauss_position = (np.random.random(1) / 5.) -3.1
    second_gauss_width = (np.random.random(1) / 10.0) + 0.05
        
    plaw_amp = (np.random.random(1) *20) - 10
    plaw_index = np.random.random(1) * (-5) - 1
    break_position = np.random.random(1) * (-4) - 1.5
    background = np.random.random(1) * (-20) 

    # return the initial guess parameters for the appropriate model
    if model == 'pow_const':
        initial_guess = [plaw_amp[0], plaw_index[0], background[0]]
            
    elif model == 'pow_const_gauss':
        initial_guess = [plaw_amp[0], plaw_index[0], background[0], gauss_amp[0],
                         gauss_position[0], gauss_width[0]]

    elif model == 'bpow_const':
        initial_guess = [plaw_amp[0], plaw_index[0], np.exp(break_position[0]), plaw_index[0], background[0]]

    elif model == 'pow_const_2gauss':
         initial_guess = [plaw_amp[0], plaw_index[0], background[0], gauss_amp[0],
                          gauss_position[0], gauss_width[0], second_gauss_amp[0],
                          second_gauss_position[0], second_gauss_width[0]]
    else:
        raise ValueError
        
    return initial_guess
    

