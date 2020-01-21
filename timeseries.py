"""
Simple time series object
"""

import numpy as np
from matplotlib import pyplot as plt


class SampleTimes:
    def __init__(self, time, label='time', units='seconds'):
        """A class holding time series sample times."""
        # ensure that the initial time is zero
        self.time = time - time[0]

        # Number of sample times
        self.nt = self.time.size

        # Average cadence
        self.dt = self.time[-1] / (self.nt - 1)

        # Information on the units for the time
        self.label = label
        self.units = units

        # Differences between consecutive sample times
        self.tdiff = self.time[1:] - self.time[0:-1]

        # Include base time for input series
        self.basetime=time[0]


class Frequencies:
    def __init__(self, frequencies, label='frequency', units='Hz'):
        self.frequencies = frequencies
        self.posindex = self.frequencies > 0
        self.positive = self.frequencies[self.posindex]
        self.label = label
        self.units = units


class PowerSpectrum:
    def __init__(self, frequencies, power, label='Fourier power'):
        self.frequencies = Frequencies(frequencies)
        self.power = power
        self.ppower = self.power[self.frequencies.posindex]
        self.label = label

        # Mean power
        self.vaughan_mean = np.mean(self.ppower)

        # Power spectrum normalized by its mean
        self.normed_by_mean = self.ppower / self.vaughan_mean

        # Standard deviation of the normalized power
        self.vaughan_std = np.std(self.normed_by_mean)

        # Normalized power expressed in units of its standard deviation
        self.Npower = self.normed_by_mean / self.vaughan_std

    def peek(self, **kwargs):
        """
        Generates a quick plot of the positive frequency part of the power
        spectrum.
        """
        plt.plot(self.frequencies.positive, self.ppower, **kwargs)


class TimeSeries:
    def __init__(self, time, data, label='data', units=None, name=None):
        self.SampleTimes = SampleTimes(time)
        if self.SampleTimes.nt != data.size:
            raise ValueError('length of sample times not the same as the data')
        self.data = data
        self.PowerSpectrum = PowerSpectrum(np.fft.fftfreq(self.SampleTimes.nt, self.SampleTimes.dt),
                                           np.abs(np.fft.fft(self.data)) ** 2)
        self.label = label
        self.units = units
        self.name = name

    def peek(self, **kwargs):
        """
        Generates a quick plot of the data
        """
        plt.plot(self.SampleTimes.time, self.data, **kwargs)
        xunits = prepend_space(bracketize(self.SampleTimes.units))
        plt.xlabel(self.SampleTimes.label + xunits)
        nsamples = ' [%i samples]' % self.SampleTimes.nt
        if self.units is not None:
            yunits = prepend_space(bracketize(self.units))
            plt.ylabel(self.label + yunits + nsamples)
        else:
            plt.ylabel(self.label + nsamples)


def prepend_left_bracket(s, bracket='(', force_replace=False,
                          test_set=('(', '{', '[')):
    """Prepend a left bracket if possible"""
    if s[0] not in test_set:
        s = bracket + s
    else:
        if force_replace:
            s[0] = bracket
    return s


def append_right_bracket(s, bracket=')', force_replace=False,
                         test_set=(')', '}', ']')):
    """Append a left bracket if possible"""
    if s[-1] not in test_set:
        s = s + bracket
    else:
        if force_replace:
            s[-1] = bracket
    return s


def bracketize(s, bracket='()', force_replace=False):
    """Add enclosing brackets if need be"""
    s = prepend_left_bracket(s, bracket=bracket[0],
                             force_replace=force_replace)
    s = append_right_bracket(s, bracket=bracket[1],
                             force_replace=force_replace)
    return s


def prepend_space(s, force_replace=False):
    s = prepend_left_bracket(s, bracket=' ', force_replace=force_replace,
                              test_set=(' '))
    return s
