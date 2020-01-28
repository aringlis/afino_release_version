import numpy as np
from afino_series import AfinoSeries


def prep_series(ts):

    # put the data series into the form (I - <I>) / <I>.
    # multiply the series by a Hanning window to aid the FFT process
    data_ave = np.mean(ts.data)
    newdata = ((ts.data - data_ave) / data_ave) * np.hanning(len(ts.data))

    ts_hann = AfinoSeries(ts.SampleTimes.time + ts.SampleTimes.basetime, newdata)
    return ts_hann
