User Guide
----------
----------

This guide will introduce you to the basics of how to perform a time series analysis of AFINO. It will also walk through the output products and how to interpret the results.

Directory structure
-------------------

By default, AFINO will save the results of an analysis in '~/afino_repository'. If this directory does not exist, it must be created. In this directory, two additional folders must be created. They are 'plots' and 'saves'.

The results can be saved in either JSON or Pickle format, via the ``use_json`` keyword. These will be placed in the 'saves' folder. AFINO will also produce a summary plot that will be placed in the 'plots' folder.

Running AFINO
-------------

The main executable function is ``analyse_series`` within the ``afino_start`` module. The example below shows a simple analysis of a sine wave in white noise::

  import numpy as np
  import afino
  from afino import afino_start

  tt = np.linspace(0,1000,1001)
  flux = np.sin(0.1 * tt) + 3*np.random.random(1001)
  afino_start.analyse_series(tt,flux, model_ids = [0,1,2,3], use_json = True) 
  
Once the analysis is complete, AFINO will generate a results file (in this example a .json file) and a summary plot. The results will identify a strong peak at P ~ 63s as expected. The default locations for these results are mentioned above, but can be changed via the ``savedir`` keyword.

Inspecting Results
------------------

We can open a results file using the JSON Python module (or the pickle module if that option was selected.::

  import json
  result = json.load(open('result = json.load(open(filename,'r'))
  
where ``filename`` corresponds to the appropriate save file. We can then inspect the contents of the file.::

  In [2]:  result.keys()
  Out[2]: dict_keys(['m0', 'm1', 'm2', 'm3'])
  
Here we can see that the results file contains a dictionary. Each key 'm0','m1' corresponds to a fit of each chosen model. For example, if only models 0 and 1 were chosen to be fit in the analysis stage, then only keys 'm0' and 'm1' would be present. In this example, models 0, 1, 2, and 3 were all attempted, so their corresponding keys are present.

The results are a nested dictionary, so each key 'm0', 'm1' etc contains its own keys.::

  In [3]: result["m0"].keys()
  Out[3]: dict_keys(['lnlike', 'model', 'BIC', 'best_fit_power_spectrum', 'frequencies', 'power', 'params', 'rchi2', 'probability', 'ID'])
  
we can see that numerous properties of the data and analysis have been saved for future inspection. We summarize these properties below.

* lnlike: The value of the best-fit log-likelihood for this model
* model: a string descriptor of the model that was fitted to the data
* BIC: the value of the Bayesian Information Criterion (BIC) for this model fit
* best_fit_power_spectrum: an array containing the values of the best fit model to the Fourier power spectrum
* frequencies: an array containing the frequencies for the Fourier power spectrum
* power: an array containing the Fourier power values
* params: A list of the best-fit model parameters. These parameters produce the best fit power spectrum.
* rchi2: The reduced chi-squared value, used to determine if this model is appropriate at all given the data
* probability: The probability value associated with the rchi2 value.
* ID: The ID number corresponding to the model used

Interpretation
--------------

To determine which of the fitted models is most appropriate for the data, we compare the Bayesian Information Criterion (BIC) values for each model. The BIC is closely related to the negative log likelihood, and is given by: -2 ln(L) + k ln(n), where L is the maximum likelihood, k is the number of free parameters, and n is the number of data points in the Fourier power spectrum.

The model that best minimizes BIC is considered the best fit to the data. For a model to be considered strongly preferred over other models, we typically require a difference in the BIC value of > 10 compared with other models. If two models have a very similar BIC then we do not have much evidence in favour of one over the other.

**For detecting oscillations**, this means we are looking for one of the models that includes an oscillation (e.g. model 1) to be strongly preferred over the models that do not contain an oscillation (e.g. model 0).

**Note**: the BIC values alone provide only a *relative* comparison between the selected models, not an absolute goodness of fit to the data. In other words, *all* of the models may fit poorly. To check this we include the additional parameter rchi2, which is a chi-squared like metric for exponentially distributed data. If the rchi2 value of a model is too large (as estimated via the associated probability value), it should be considered inappropriate for the data. Suggested cutoff values for the probability are p < 0.05 or p < 0.01.


