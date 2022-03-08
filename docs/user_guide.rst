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





