import numpy as np

def rebin_prototype(orig_array,ratio):
    old_shape=np.float(orig_array.shape[0])    
    new_shape=old_shape/ratio
    #test if old_shape divided evenly
    if round(new_shape) != new_shape:
        print('new array not an integer division of the old array. End of new array will be truncated.')
    else:
        print('OK: new array is integer division of old array.')
    
    new_array=[]
    for n in range(0,np.int(old_shape),ratio):
        if n+ratio > old_shape:
            break
        else:
            tmp=np.mean(orig_array[n:n+ratio])
            new_array.append(tmp)

    result=np.array(new_array)
    return result
