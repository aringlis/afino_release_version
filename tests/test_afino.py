# tests for AFINO
import numpy as np
from afino import afino_series


def test_afinoseries():
    tt = np.linspace(0,100,101)
    flux = np.random.random(101)
    ts = afino_series.AfinoSeries(tt,flux)

    assert ts.SampleTimes.nt == 101
    assert ts.SampleTimes.dt == 1.0
    assert len(ts.data) == 101
    assert len(ts.PowerSpectrum.frequencies.positive) == 50
    assert len(ts.PowerSpectrum.ppower) == 50

def test_prep_series():
    tt = np.linspace(0,100,101)
    flux = np.random.random(101)
    ts = afino_series.AfinoSeries(tt,flux)
    ts_prepped = afino_series.prep_series(ts)

    assert ts_prepped.data[0] == 0.0
    assert ts_prepped.data[-1] == 0.0

def test_nothing():
    pass



    

