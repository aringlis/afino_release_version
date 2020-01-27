import pickle
import os
import numpy as np
from afino_main_analysis3 import main_analysis


def model_comparison(ts,description=None,generic=False,low_frequency_cutoff=None, overwrite_gauss_bounds = None):


    m0 = main_analysis(ts, model='single_power_law_with_constant',low_frequency_cutoff=low_frequency_cutoff)
  
    m1 = main_analysis(ts, model='splwc_AddNormalBump2',low_frequency_cutoff=low_frequency_cutoff, overwrite_gauss_bounds = overwrite_gauss_bounds)

    m2 = main_analysis(ts, model='broken_power_law_with_constant',low_frequency_cutoff=low_frequency_cutoff)

    #parse everything into one dictionary

    dBIC_0v1 = m0['BIC'] - m1['BIC']
    dBIC_2v1 = m2['BIC'] - m1['BIC']
    dBIC_0v2 = m0['BIC'] - m2['BIC'] 
    analysis_summary = {}
    analysis_summary['m0'] = m0
    analysis_summary['m1'] = m1
    analysis_summary['m2'] = m2
    analysis_summary['dBIC'] = dBIC_0v1
    analysis_summary['dBIC_2v1'] = dBIC_2v1
    analysis_summary['dBIC_0v2'] = dBIC_0v2

    #print some info for the user
    print('Analysis summary info:')
    print('-----------------------------')
    print('Lnlike plaw: ' + str(m0['lnlike']))
    print('Lnlike gauss: ' + str(m1['lnlike']))
    print('Lnlike bpow: ' + str(m2['lnlike']))
    print(' ')
    print('BIC plaw: ' + str(m0['BIC']))
    print('BIC gauss: ' + str(m1['BIC']))
    print('BIC bpow: ' + str(m2['BIC']))
    print(' ')
    print('rchi2 plaw: ' + str(m0['rchi2']))
    print('prob. plaw: ' + str(m0['probability']))
    print(' ')
    print('rchi2 gauss: ' + str(m1['rchi2']))
    print('prob. gauss: ' + str(m1['probability']))
    print(' ')
    print('rchi2 bpow: ' + str(m2['rchi2']))
    print('prob. bpow: ' + str(m2['probability']))
    print(' ')
    print('dBIC M0 vs M1: ' + str(dBIC_0v1))
    print('dBIC M2 vs M1: ' + str(dBIC_2v1))
    print('dBIC M0 vs M2: ' + str(dBIC_0v2))
    print('-----------------------------')
    print(' ')
    print(' ')
    

    #now have a dictionary of results for each model in m0 and m1
    #save this to a combined results file

    if generic:
        fname = os.path.join(os.path.expanduser('~/afino_repository/generic/'),'afino_summary_' + description + '.pickle')
    else:
        fname = os.path.join(os.path.expanduser('~/afino_repository/afino_analysis_summaries/'),'afino_summary_' + description + '.pickle')
    pickle.dump(analysis_summary,open(fname,'wb'))
    analysis_summary['analysis_filename'] = fname

    return analysis_summary
    
   
  

def model_comparison_mms(ts,root='~/proba2_results',description='tmp',generic=False,low_frequency_cutoff=None, overwrite_gauss_bounds = None, overwrite_extra_gauss_bounds = None):


    m0 = main_analysis(ts, model='single_power_law_with_constant',low_frequency_cutoff=low_frequency_cutoff)
  
    m1 = main_analysis(ts, model='splwc_AddNormalBump2',low_frequency_cutoff=low_frequency_cutoff, overwrite_gauss_bounds = overwrite_gauss_bounds)

    m2 = main_analysis(ts, model='broken_power_law_with_constant',low_frequency_cutoff=low_frequency_cutoff)

    m3 = main_analysis(ts, model='splwc_AddNormalBump2_plus_extra_bump',low_frequency_cutoff=low_frequency_cutoff, overwrite_extra_gauss_bounds = overwrite_extra_gauss_bounds)

    m4 = main_analysis(ts, model='splwc_AddNormalBump2',low_frequency_cutoff = low_frequency_cutoff, overwrite_gauss_bounds = overwrite_extra_gauss_bounds[0:3] + overwrite_extra_gauss_bounds[6:9])

    #parse everything into one dictionar

    dBIC_0v1 = m0['BIC'] - m1['BIC']
    dBIC_2v1 = m2['BIC'] - m1['BIC']
    dBIC_3v1 = m3['BIC'] - m1['BIC']
    dBIC_4v1 = m4['BIC'] - m1['BIC']
    dBIC_0v2 = m0['BIC'] - m2['BIC']
    dBIC_3v2 = m3['BIC'] - m2['BIC']
    dBIC_4v2 = m4['BIC'] - m2['BIC']
    dBIC_0v3 = m0['BIC'] - m3['BIC']
    dBIC_4v3 = m4['BIC'] - m3['BIC']
    dBIC_0v4 = m0['BIC'] - m4['BIC']
    analysis_summary = {}
    analysis_summary['m0'] = m0
    analysis_summary['m1'] = m1
    analysis_summary['m2'] = m2
    analysis_summary['m3'] = m3
    analysis_summary['m4'] = m4
    analysis_summary['dBIC'] = dBIC_0v1
    analysis_summary['dBIC_2v1'] = dBIC_2v1
    analysis_summary['dBIC_0v2'] = dBIC_0v2
    analysis_summary['dBIC_3v1'] = dBIC_3v1
    analysis_summary['dBIC_3v2'] = dBIC_3v2
    analysis_summary['dBIC_0v3'] = dBIC_0v3
    analysis_summary['dBIC_4v1'] = dBIC_4v1
    analysis_summary['dBIC_4v2'] = dBIC_4v2
    analysis_summary['dBIC_4v3'] = dBIC_4v3
    analysis_summary['dBIC_0v4'] = dBIC_0v4

    #print some info for the user
    print('Analysis summary info:')
    print('-----------------------------')
    print('Lnlike plaw: ' + str(m0['lnlike']))
    print('Lnlike gauss: ' + str(m1['lnlike']))
    print('Lnlike bpow: ' + str(m2['lnlike']))
    print(' ')
    print('BIC plaw: ' + str(m0['BIC']))
    print('BIC gauss: ' + str(m1['BIC']))
    print('BIC bpow: ' + str(m2['BIC']))
    print('BIC gauss_plus_extra: ' + str(m3['BIC']))
    print(' ')
    print('rchi2 plaw: ' + str(m0['rchi2']))
    print('prob. plaw: ' + str(m0['probability']))
    print(' ')
    print('rchi2 gauss: ' + str(m1['rchi2']))
    print('prob. gauss: ' + str(m1['probability']))
    print(' ')
    print('rchi2 bpow: ' + str(m2['rchi2']))
    print('prob. bpow: ' + str(m2['probability']))
    print(' ')
    print('dBIC M0 vs M1: ' + str(dBIC_0v1))
    print('dBIC M2 vs M1: ' + str(dBIC_2v1))
    print('dBIC M0 vs M2: ' + str(dBIC_0v2))
    print('-----------------------------')
    print(' ')
    print(' ')
    

    #now have a dictionary of results for each model in m0 and m1
    #save this to a combined results file

    if generic:
        fname = os.path.join(os.path.expanduser('~/afino_repository/generic/'),'afino_summary' + description + '.pickle')
    else:
        fname = os.path.join(os.path.expanduser('~/afino_repository/afino_analysis_summaries/'),'afino_summary' + description + '.pickle')
    pickle.dump(analysis_summary,open(fname,'wb'))
    analysis_summary['analysis_filename'] = fname

    return analysis_summary
    
