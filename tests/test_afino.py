# tests for AFINO
import numpy as np
from numpy.testing import assert_almost_equal
from afino import afino_series
from afino import afino_utils
from afino import afino_spectral_models
from afino.afino_main_analysis3 import main_analysis

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

def test_main_analysis():
    tt = np.linspace(0,100,101)
    flux = np.random.random(101)
    ts = afino_series.AfinoSeries(tt,flux)
    ts_prepped = afino_series.prep_series(ts)
    result = main_analysis(ts_prepped, model='pow_const')

    assert type(result) == dict
    assert type(result['lnlike']) == np.float64
    assert type(result['BIC']) == np.float64
    assert result['model'] == 'pow_const'
    assert type(result['frequencies']) == np.ndarray
    assert type(result['power']) == np.ndarray
    assert type(result['best_fit_power_spectrum']) == np.ndarray
    
# various tests for the model functions

def test_pow():
    freqs = np.linspace(0.001,0.5)
    power = afino_spectral_models.pow([1.0,1.0],freqs)
    assert_almost_equal(power[0], 2.718, decimal = 3)

def test_pow_const():
    freqs = np.linspace(0.001,5.0)
    power = afino_spectral_models.pow_const([1.0,1.0,0.0],freqs)
    assert_almost_equal(power[-1], 1.00, decimal = 2)
    
    



    

