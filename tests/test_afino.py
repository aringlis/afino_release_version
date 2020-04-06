# tests for AFINO
import numpy as np
from numpy.testing import assert_almost_equal
from afino import afino_series
from afino import afino_utils
from afino import afino_spectral_models


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

def test_model_id_to_string():
    model_string = afino_utils.model_string_from_id(1)
    assert model_string == 'pow_const_gauss'
    model_string = afino_utils.model_string_from_id(3)
    assert model_string == 'pow_const_2gauss'

def test_nothing():
    pass


# various tests for the model functions

def test_pow():
    freqs = np.linspace(0.001,0.5)
    power = afino_spectral_models.pow([1.0,1.0],f)
    assert_almost_equal(power[0], 2.718, decimal = 3)

    



    

