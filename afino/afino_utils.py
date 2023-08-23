"""
This module provides a number of convenience functions for use throughout
the AFINO codebase.
"""


import pickle
import json
import os
import copy
import numpy as np
import itertools

class NumpyEncoder(json.JSONEncoder):
    """This encoder converts numpy arrays to lists so that they
       can be saved in JSON format"""
    
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


    
def model_string_from_id(id):
    """This function converts an AFINO model ID integer to a string descriptor"""

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


def relative_bics(saveresult):
    """
    This convenience function calculates all the relative BIC comparisons
    from an AFINO save result.

    Parameters
    ----------
    saveresult : dict
        An AFINO save result dictionary. Can be restored from a previous save file

    Returns
    -------
    dbic_values : dict
        A dictionary of relative BIC value comparisons, i.e. BICa - BICb for every a,b pair
    
    """
    
    bics = []
    ids = []   
    dbic_values = {}

    for key in saveresult:
        bics.append(saveresult[key]['BIC'])
        ids.append(int(saveresult[key]['ID']))

    # extract all the ID and BIC pairs (includes reverse pairings)
    combo_list = list(itertools.permutations(ids))
    combo_bic_list = list(itertools.permutations(bics))

    # for each ID/BIC pair a,b, put BICa - BICb in a dictionary
    for i, combo in enumerate(combo_list):
        dbic_values['dBIC_' + str(combo[0]) + '_minus_' + str(combo[1])] = combo_bic_list[i][0] - combo_bic_list[i][1]
    
    return dbic_values
    

def save_afino_results(results, use_json = False, description = None, nosave = False):
    """
    This function saves the results of an AFINO analysis run to either a JSON file
    or a Pickle file.

    Parameters
    ----------

    results : list
        An AFINO results list, obtained via the model_comparison function
    use_json : bool, optional
        If True, use JSON format to save the analysis results. If False, use pickle format
    description : string, optional
        An optional descriptor for the result that is incorporated into the filename

    """
    

    # ensure the directory to save plots exists, create it if not.
    os.makedirs(os.path.expanduser('~/afino_repository/saves/'),exist_ok=True)
    
    analysis_summary = {}

    # This may not be very Pythonic
    for i, r in enumerate(results):
        id = r['ID']
        analysis_summary['m'+str(id)] = r

    # print some info for the user
    print(' ')
    print('Analysis summary info:')
    print('-----------------------------')
    print(' ')
    for r in results:
        print('Lnlike m' + str(r['ID']) + ' (' + r['model'] + '): ' + str(r['lnlike']))
    print(' ')
    for r in results:
        print('BIC m' + str(r['ID']) + ' (' + r['model'] + '): ' + str(r['BIC']))
    print(' ')
    for r in results:
        print('rchi2 m' + str(r['ID']) + ' (' + r['model'] + '): ' + str(r['rchi2']))
        print('prob. m' + str(r['ID']) + ' (' + r['model'] + '): ' + str(r['probability']))
    print(' ')
    print('-----------------------------')
    print(' ')
    print(' ')
        

    # save all the results to a JSON or pickle file

    if not nosave:
    
        if use_json:
            fname = os.path.join(os.path.expanduser('~/afino_repository/saves/'),'afino_summary_data_' + description + '.json')
            json.dump(analysis_summary,open(fname,'w'), cls = NumpyEncoder)
        else:
            fname = os.path.join(os.path.expanduser('~/afino_repository/saves/'),'afino_summary_data_' + description + '.pickle')
            pickle.dump(analysis_summary,open(fname,'wb'))
        

    return analysis_summary

        

def restore_json_save_file(fname):
    """This function restores an AFINO JSON save file and converts
       certain dictionary entries back into their original ndarray form.
       The output should be identical to that from the pickle save files"""
    
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
    

    
