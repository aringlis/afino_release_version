

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
        
    
