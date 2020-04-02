
import afino_start
from astropy.io import fits

result = fits.open('test_data/flare566801.fits')
hdu1 = result[1]
flux = hdu1.data['flux']
tt = hdu1.data['time']
tt2 = tt*60
afino_start.analyse_series(tt2,flux,low_frequency_cutoff=10,
                                      overwrite_gauss_bounds = [(-10.0,10.0),(-1.0,6.0),
                                                                (-20.0,10.0),(-16.0,5.0),(-9.0,-5.0),(0.05,0.25)]) 
