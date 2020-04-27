#
# Fit an arbitrary function with a model using maximum likelihood,
# assuming exponential distributions.
#
import numpy as np
import scipy.optimize as op
from scipy.stats import gamma
from scipy.special import gammaincc, gammainccinv


#
# Log likelihood function.  In this case we want the product of exponential
# distributions.
#
def lnlike(variables, x, y, model_function):
    """
    Log likelihood of the data given a model.  Assumes that the data is
    exponentially distributed.  Can be used to fit Fourier power spectra.
    :param variables: array like, variables used by model_function
    :param x: the independent variable (most often normalized frequency)
    :param y: the dependent variable (observed power spectrum)
    :param model_function: the model that we are using to fit the power
    spectrum
    :return: the log likelihood of the data given the model.
    """
    model = model_function(variables, x)
    return -np.sum(np.log(model)) - np.sum(y / model)


#
# Fit the input model to the data.
#
def go_gauss(freqs, data, model_function, initial_guess, method, overwrite_gauss_bounds = None):
    nll = lambda *args: -lnlike(*args)
    args = (freqs, data, model_function)
    if overwrite_gauss_bounds:
        return op.minimize(nll, initial_guess, args=args, method=method, bounds = overwrite_gauss_bounds)
    else:
        return op.minimize(nll, initial_guess, args=args, method=method, bounds = [(-10.0,10.0),(-1.0,6.0),(-20.0,10.0),(-16.0,5.0),(-5.7,-1.5),(0.05,0.25)])

def go_plaw(freqs, data, model_function, initial_guess, method):
    nll = lambda *args: -lnlike(*args)
    args = (freqs, data, model_function)
    return op.minimize(nll, initial_guess, args=args, method=method, bounds = [(-10.0,10.0),(-1.0,6.0),(-20.0,10.0)])


def go_bpow(freqs, data, model_function, initial_guess, method):
    nll = lambda *args: -lnlike(*args)
    args = (freqs, data, model_function)
    return op.minimize(nll, initial_guess, args=args, method=method, bounds = [(None,None),(1.0,9.0),(0.0033,0.25),(1.0, 9.0),(None,None)])

def go_gauss_plus_extra_bump(freqs, data, model_function, initial_guess, method, overwrite_extra_gauss_bounds=None):
    nll = lambda *args: -lnlike(*args)
    args = (freqs, data, model_function)
    if overwrite_extra_gauss_bounds:
        return op.minimize(nll, initial_guess, args=args, method=method, bounds = overwrite_extra_gauss_bounds)
    else:
        return op.minimize(nll, initial_guess, args=args, method=method, bounds = [(-10.0,10.0),(-1.0,6.0),(-20.0,10.0),(-16.0,5.0),(-5.7,-1.5),(0.05,0.25), (-16.0,5.0),(-3.1,-2.9),(0.05,0.12)])
#bounds = [(-10.0,10.0),(-1.0,9.0),(0.01,None),(-1.0, 9.0),(-20.0,10.0)])

#
# The code below refers to equations in Nita et al (2014), ApJ, 789, 152
#
#
# Sample to Model Ratio (SMR) estimator
#
def rhoj(Sj, shatj):
    """
    Sample to Model Ratio (SMR) estimator (Eq. 5)
    :param Sj: random variables (data)
    :param shatj: best estimate of the model
    :return: SMR estimator
    """
    return Sj / shatj


#
# Goodness-of-fit estimator
#
def rchi2(m, nu, rhoj):
    """
    Goodness-of-fit estimator (Eq. 16)

    Parameters
    ----------

    m
        number of spectra considered
    nu
        degrees of freedom
    rhoj
        sample to model ratio estimator

    Returns
    -------
    float
        the sample-to-model ratio, a goodness of fit estimator
    """
    return (m / (1.0 * nu)) * np.sum((1.0 - rhoj) ** 2)


#
# PDF of the goodness-of-fit estimator (Eq. 17)
#
def rchi2distrib(m, nu):
    """
    The distribution of rchi2 may be approximated by the analytical expression
    below.  Comparing Eq. (2) with the implementation in scipy stats we find
    the following equivalencies:
    k (Nita parameter) = a (scipy stats parameter)
    theta (Nita parameter) = 1 / lambda (scipy stats value)
    :param m: number of spectra considered
    :param nu: degrees of freedom
    :return: a frozen scipy stats function that represents the distribution
    of the data
    """
    # Calculate the gamma function parameter values as expressed in Eq. 17
    k = (nu / 2.0) / (1.0 + 3.0 / m)
    scale = 1.0 / k
    #
    return gamma(k, scale=scale, loc=0.0)


#
# Probability of getting this value of reduced chi-squared or larger (Eq. 18)
#
def prob_this_rchi2_or_larger(rchi2, m, nu):
    """
    :param rchi2: reduced chi-squared value
    :param m: number of spectra considered
    :param nu:  degrees of freedom
    :return:
    """
    a = (nu / 2.0) * np.float64(m) / (3.0 + np.float64(m))
    return gammaincc(a, a * rchi2)


#
# What value of rchi2 gives rise to a given probability level?
#
def rchi2_given_prob(p, m, nu):
    a = (nu / 2.0) * np.float64(m) / (3.0 + np.float64(m))
    return gammainccinv(a, p) / a


def AIC(k, variables, freqs, data, model_function):
    return 2 * k - 2 * lnlike(variables, freqs, data, model_function)


def BIC(k, variables, freqs, data, model_function, n):
    return -2 * lnlike(variables, freqs, data, model_function) + k * np.log(n)
