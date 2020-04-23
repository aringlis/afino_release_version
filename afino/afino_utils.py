import pickle
import json
import os
import copy
import numpy as np

class NumpyEncoder(json.JSONEncoder):
    '''This encoder converts numpy arrays to lists so that they
       can be saved in JSON format'''
    
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


    
def model_string_from_id(id):

    allowed_models = {
        0 : 'pow_const',
        1 : 'pow_const_gauss',
        2 : 'bpow_const',
        3 : 'pow_const_2gauss'
    }

    model_string = allowed_models.get(id)
    if not model_string:
        raise ValueError('Invalid model string')
    
    return model_string
        
    

def save_afino_results(results, use_json = False, description = None):

    #ensure the directory to save plots exists, create it if not.
    os.makedirs(os.path.expanduser('~/afino_repository/saves/'),exist_ok=True)
    
    analysis_summary = {}

    #almost certainly not pythonic!
    for i, r in enumerate(results):
        analysis_summary['m'+str(i)] = r

    m0 = results[0]
    m1 = results[1]
    m2 = results[2]
    #add extras to dictionary - need to change this
    dBIC_0v1 = results[0]['BIC'] - results[1]['BIC']
    dBIC_2v1 = results[2]['BIC'] - results[1]['BIC']
    dBIC_0v2 = results[0]['BIC'] - results[2]['BIC'] 

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

    #save all the results to a JSON or pickle file

    if use_json:
        fname = os.path.join(os.path.expanduser('~/afino_repository/saves/'),'afino_summary_data_' + description + '.json')
        json.dump(analysis_summary,open(fname,'w'), cls = NumpyEncoder)
    else:
        fname = os.path.join(os.path.expanduser('~/afino_repository/saves/'),'afino_summary_data_' + description + '.pickle')
        pickle.dump(analysis_summary,open(fname,'wb'))
        

    return analysis_summary


def restore_json_save_file(fname):
    '''This function restores an AFINO JSON save file and converts
       certain dictionary entries back into their original ndarray form.
       The output should be identical to that from the pickle save files'''
    
    result = json.load(open(fname,'r'))
    converted_result = copy.deepcopy(result)

    # iterate through the dictionary keys and convert the best fit power spectrum,
    # frequencies, and power keys back into numpy arrays
    keys_to_convert = ['best_fit_power_spectrum', 'frequencies', 'power']
    for key in result.keys():
        if isinstance(result[key],dict):
            for subkey in result[key].keys():
                if subkey in keys_to_convert:
                    converted_result[key][subkey] = np.asarray(result[key][subkey])
                    
    return converted_result
    

    
