"""
Power Spectrum Models
"""
from copy import deepcopy
import numpy as np
import astropy.units as u
#import lnlike_model_fit
#import pstools
import matplotlib.pyplot as plt


#
# Normalize the frequency
#
def fnorm(f, normalization):
    """ Normalize the frequency spectrum."""
    return f / normalization


# ----------------------------------------------------------------------------
# constant
#
def constant(a):
    """The power spectrum is a constant across all frequencies
    
    Parameters
    ----------
    a : float
        the natural logarithm of the power
    """
    return np.exp(a)


# ----------------------------------------------------------------------------
# Power law
#
def pow(a, f):
    """Simple power law.  This model assumes that the power
    spectrum is made up of a power law at all frequencies.
    
    Parameters
    ----------
    a : ndarray[2]
        a[0] : the natural logarithm of the normalization constant
        a[1] : the power law index
    f : ndarray
        frequencies
    """
    return np.exp(a[0]) * ((fnorm(f, f[0]) ** (-a[1])))# (f ** (-a[1]))


# ----------------------------------------------------------------------------
# Broken Power law
#
def bpow(a, f):
    """Broken power law.  This model assumes that there is a break in the power
    spectrum at some given frequency.
   
    Parameters
    ----------
    a : ndarray[3]
        a[0] : the natural logarithm of the normalization constant
        a[1] : the power law index at frequencies lower than the break frequency
        a[2] : break frequency
        a[3] : the power law index at frequencies higher than the break frequency
    f : ndarray
        frequencies
    """
    power = np.zeros_like(f)
    less_than_break = f < a[2]
    above_break = f >= a[2]
    power[less_than_break] = np.exp(a[0]) * f[less_than_break] ** (-a[1])
    power[above_break] = np.exp(a[0]) * (a[2]**(-a[1]+a[3])) * f[above_break] ** (-a[3])
    return power


# ----------------------------------------------------------------------------
# Sum of Pulses
#
def sum_of_pulses(a, f):
    """Sum of pulses.  This model is based on Aschwanden "Self Organized
    Criticality in Astrophysics", Eq. 4.8.23.  Simulations implementing this
    equation come up with a shape that is modeled below.
    
    Parameters
    ----------
    a : ndarray[3]
        a[0] : the natural logarithm of the normalization constant
        a[1] : the scale frequency
        a[2] : the power law index
    f : ndarray
        frequencies
    """
    return np.exp(a[0])/(1.0 + (f/a[1]) ** a[2])


# ----------------------------------------------------------------------------
# Sum of Pulses Plus Constant
#
def sum_of_pulses_with_constant(a, f):
    """Sum of pulses plus constant.  This model is based on Aschwanden "Self
    Organized Criticality in Astrophysics", Eq. 4.8.23, with a constant
    background to model detector noise.
    
    Parameters
    ----------
    a : ndarray[3]
        a[0] : the natural logarithm of the normalization constant
        a[1] : the scale frequency
        a[2] : the power law index
        a[3] : natural logarithm of the background constant
    f : ndarray
        frequencies
    """
    return np.exp(a[0])/(1.0 + (f/a[1]) ** a[2]) + np.exp(a[3])


# ----------------------------------------------------------------------------
# Broken Power law with Constant
#
def bpow_const(a, f):
    """Broken power law with constant.  This model assumes that there is a
    break in the power spectrum at some given frequency.  At high
    frequencies the power spectrum is dominated by the constant background.
    
    Parameters
    ----------
    a : ndarray(5)
        a[0] : the natural logarithm of the normalization constant
        a[1] : the power law index at frequencies lower than the break frequency
        a[2] : break frequency
        a[3] : the power law index at frequencies higher than the break frequency
        a[4] : the natural logarithm of the constant background
    f : ndarray
        frequencies
    """
    return bpow(a[0:4], f) + constant(a[4])


def broken_power_law_with_constant_plus_extra_bump(a, f):

    return bpow(a[0:4], f) + constant(a[4]) + NormalBump2(np.log(f), a[5:8])
    

# ----------------------------------------------------------------------------
# Broken Power law with Constant
#
def broken_power_law_with_constant_with_lognormal(a, f):
    """Broken power law with constant with lognormal.  This model assumes that
    there is a break in the power spectrum at some given frequency.  At high
    frequencies the power spectrum is dominated by the constant background.  At
    some particular frequency there is a lognormal (narrowband distribution)
    
    Parameters
    ----------
    a : ndarray(5)
        a[0] : the natural logarithm of the normalization constant
        a[1] : the power law index at frequencies lower than the break frequency
        a[2] : break frequency
        a[3] : the power law index at frequencies higher than the break frequency
        a[4] : the natural logarithm of the constant background
    f : ndarray
        frequencies
    """
    return bpow(a[0:4], f) + constant(a[4])


# ----------------------------------------------------------------------------
# Power law with constant
#
def pow_const(a, f):
    """Power law with a constant.  This model assumes that the power
    spectrum is made up of a power law and a constant background.  At high
    frequencies the power spectrum is dominated by the constant background.
    
    Parameters
    ----------
    a : ndarray[2]
        a[0] : the natural logarithm of the normalization constant
        a[1] : the power law index
        a[2] : the natural logarithm of the constant background
    f : ndarray
        frequencies
    """
    return pow(a[0:2], f) + constant(a[2])


# Lognormal
#
def lognormal(a, f):
    """
    A lognormal distribution
    
    Parameters
    ----------
    a : ndarray(3)
        a[0] : the natural logarithm of the Gaussian amplitude
        a[1] : the natural logarithm of the center of the Gaussian
        a[2] : the width of the Gaussian in units of natural logarithm of the
               frequency
    f : ndarray
        frequencies
    """
    onent = (np.log(f) - a[1]) / a[2]
    amp = np.exp(a[0])
    return amp * np.exp(-0.5 * onent ** 2)


def lognormal_CF(f, a, b, c):
    return lognormal([a, b, c], f)


def power_law_with_constant_with_lognormal(a, f):
    """
    Power law with constant and a lognormal.
    
    Parameters
    ----------
    a : ndarray[6]
        a[0] : the natural logarithm of the normalization constant
        a[1] : the power law index
        a[2] : the natural logarithm of the constant background
        a[3] : the natural logarithm of the Gaussian amplitude
        a[4] : the natural logarithm of the center of the Gaussian
        a[5] : the width of the Gaussian in units of natural logarithm of the
               frequency
    f : ndarray
        frequencies
    """
    return pow_const(a[0:3], f) + lognormal(a[3:6], f)







def NormalBump2(x, a):
    z = (x - a[1]) / a[2]
    amplitude = np.exp(a[0])
    norm = 1.0 / (np.sqrt(2 * np.pi * a[2] ** 2))
    return amplitude * norm * np.exp(-0.5 * z ** 2)

def pow_const_gauss(a, f):
    """Simple power law with a constant, plus a Gaussian shaped bump.
    This model assumes that the powe spectrum is made up of a power law and a
    constant background.  At high frequencies the power spectrum is dominated
    by the constant background.

    Parameters
    ----------
    f : ndarray
        frequencies

    a : ndarray[2]
        a[0] : the natural logarithm of the normalization constant
        a[1] : the power law index
        a[2] : the natural logarithm of the constant background
    """
    return pow_const(a[0:3], f) + NormalBump2(np.log(f), a[3:6])


def pow_const_2gauss(a,f):
    """Simple power law with a constant, plus a Gaussian bump, and a second Gaussian bump.
    The second bump is implemented to account for a persistent signal in the data, such as a spin period.

    Parameters
    ----------
    f : ndarray
        frequencies

    a : ndarray[8]
        a[0] : the natural logarithm of the normalization constant
        a[1] : the power law index
        a[2] : the natural logarithm of the constant background
        a[3] : The amplitude of the first bump
        a[4] : the frequency location of the first bump
        a[5] : the width of the first bump
        a[6] : The amplitude of the second bump
        a[7] : the frequency location of the second bump
        a[8] : the width of the second bump
    """
        
    return pow_const(a[0:3], f) + NormalBump2(np.log(f), a[3:6]) + NormalBump2(np.log(f), a[6:9])

    
