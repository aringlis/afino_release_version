# tests for AFINO
import numpy as np
from afino import afino_series


def test_prep_series():
    tt = np.linspace(0,100,101)
    flux = np.random.random(101)
    ts = afino_series.AfinoSeries(tt,flux)
    ts_prepped = afino_series.prep_series(ts)

    assert ts_prepped.data[0] == 0.0
    assert ts_prepped.data[-1] == 0.0

def test_nothing():
    pass



    

