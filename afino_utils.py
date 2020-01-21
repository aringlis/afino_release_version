#!/Users/ainglis/anaconda/bin/python



import os
import pickle
import sys
import csv
import glob
import matplotlib.pyplot as plt
sys.path.append('/Users/Inglis/physics/rednoise/py/tools/')
sys.path.append('/Users/Inglis/physics/rednoise/py/testing/')
from afino_main_analysis import main_analysis
import rnspectralmodels3
from afino_model_comparison import model_comparison
from load_event_data import prep_series
from timeseries import TimeSeries
import shutil
import datetime
import urllib.parse
import numpy as np

import sunpy
from sunpy import lightcurve
from sunpy.net import hek
from sunpy.time import parse_time,TimeRange
from goes_event_list import get_goes_event_list
from sunpy.instr import rhessi, lyra

from astropy.io import fits
import afino_gbm_utils

#define AFINO directory paths as global variables
ANALYSIS_RECORD_DIRECTORY = os.path.expanduser('~/afino_repository/analysis_records')
PERMANENT_ANALYSIS_RESULT_DIRECTORY = os.path.expanduser('~/afino_repository/mcmc_analysis_results')
MODEL_COMPARISON_DIRECTORY=os.path.expanduser('~/afino_repository/model_comparisons')
ANALYSIS_SUMMARY_PLOTS_DIRECTORY = os.path.expanduser('~/afino_repository/analysis_summary_plots')
AFINO_DATA_DIRECTORY = os.path.expanduser('~/afino_repository/data')
AFINO_EVENT_LIST_DIRECTORY = os.path.expanduser('~/afino_repository/event_lists')
POSTERIOR_PREDICTIVE_DIRECTORY = os.path.expanduser('~/afino_repository/posterior_predictive_results')

def search_for_events_on_date(date=None):
    '''This function searches the GOES event list for solar flares occuring on the given date.
    If date is not specified, then yesterday's date is used for the search.
    The default threshold for selecting an event is > GOES-class M1. The relevant GBM, LYRA and RHESSI data are downloaded.
    This function also creates a pickle file containing a list of the events of interest.

    Data are downloaded to ~/afino_repository/data/
    Event lists are downloaded to ~/afino_repository/event_lists.'''


    
    now=datetime.datetime.utcnow()
    date_today=now.date()
    dt=datetime.timedelta(days=1)
    date_yesterday=date_today-dt
    if not date:
        #if date was not passed in then search for events from yesterday
        tran=TimeRange(date_yesterday.isoformat(),date_today.isoformat())
    else:
        #if date was passed in then make a time range spanning that day
        date2=parse_time(date)
        tran=TimeRange(date2.date().isoformat(),(date2+datetime.timedelta(1)).date().isoformat())

    #get any events that occured in the previous day above GOES-class M1
    #events is a list of dictionaries, each one containing information on a flare
    events=get_goes_event_list(tran,goes_class_filter='M1')
    n=len(events)

    if n == 0:
        print('No events found in range, no data will be downloaded')
    
    #print some info about what the search found
    print(' ')
    print(str(n) +' events matching the criteria occured between '+tran.start.date().isoformat() + ' and '+tran.end.date().isoformat() +':')
    print(' ')
    print('----------------------------------------------------')
    for i, e in enumerate(events):
        print('Event '+str(i) +':  '+ e['start_time'].isoformat() + ' - ' + e['end_time'].isoformat() + '   GOES class: ' + e['goes_class'])
    print('----------------------------------------------------')
    print(' ')



    #next step is to download the data for the events and instruments in question.
    #tools already exist to download RHESSI obs_summ data in SunPy

    dir=AFINO_DATA_DIRECTORY #os.path.expanduser('~/afino_repository/data/')
    if not os.path.exists(dir):
        os.makedirs(dir)

    #for ev in events:
    if n > 0:

        #-------------------------------------------------------------------------------
    
        #download the RHESSI data using SunPy tools
       # print 'Getting RHESSI data file(s)...'
       # print '----------------------------'
        ev_tran=tran #TimeRange(ev['start_time'],ev['end_time'])
       # rhessi_file_name=rhessi.get_obssum_filename(ev_tran)
       # rhessi_file=rhessi.get_obssumm_file(ev_tran)
       # destination_rhessi_file_name=str.split(urlparse.urlsplit(rhessi_file_name).path,'/')[-1]
       # shutil.copy(rhessi_file[0],os.path.join(dir,destination_rhessi_file_name))
       # print 'Moved RHESSI file from '+ rhessi_file[0] +' to '+os.path.join(dir,rhessi_file_name)
        #rhessi_data_dict=rhessi.parse_obssumm_file(rhessi_file[0])

        #-------------------------------------------------------------------------------
    
        #also want to download the lyra data. Use a lyra lightcurve to do this.
       # print ' '
       # print 'Getting LYRA data file(s)...'
       # print '----------------------------'
       # l=lightcurve.LYRALightCurve.create(ev_tran.start.date().isoformat(), overwrite=True)
    
        #data file is downloaded to ~/sunpy/data/. Copy to a better location
        #do some string manipulation to get the file name
        datestr=str.split(ev_tran.start.date().isoformat(),'-')
       # fname_datestring=datestr[0]+datestr[1]+datestr[2]
    
       # fname='lyra_'+fname_datestring+'-000000_lev2_std.fits'
        #copy the file to dir
       # shutil.copy(os.path.join(os.path.expanduser('~/sunpy/data/'),fname), dir)
       # print 'Moved LYRA file from ~/sunpy/data/ to ' + dir

        #-------------------------------------------------------------------------------
    
        #also want to download Fermi/GBM data. Use the Fermi/GBM lightcurve to do this.
       # print ' '
       # print 'Getting Fermi/GBM data file(s)...'
       # print '--------------------------------'
       # det='n1'
       # gbm_l=lightcurve.GBMSummaryLightCurve.create(ev_tran.start.date().isoformat(),detector=det)

        #data file is downloaded to ~/sunpy/data./ Copy to a better location
        #do some string manipulation to get the file name
       # fermi_datestring=datestr[0][2:4] + datestr[1] + datestr[2]
       # fname='glg_cspec_' + det + '_' + fermi_datestring + '_v00.pha'
        #copy the file to dir
       # shutil.copy(os.path.join(os.path.expanduser('~/sunpy/data/'),fname), dir)
       # print 'Moved GBM file from ~/sunpy/data to ' + dir 

        #--------------------------------------------------------------------------------
    
        #also want to download GOES data. Use the GOES lightcurve to do this.
        print(' ')
        print('Getting GOES data file(s)...')
        print('--------------------------------')
        goes_l = lightcurve.GOESLightCurve.create(ev_tran)
        goes_satnum = goes_l.meta['TELESCOP'].split()[1]
        goes_fname = 'go' + goes_satnum + datestr[0] + datestr[1] + datestr[2] + '.fits'
        #copy the file to dir
        shutil.copy(os.path.join(os.path.expanduser('~/sunpy/data/'),goes_fname), dir)
        print('Moved GOES file from ~/sunpy/data to ' + dir) 
        
        
    datestr = str.split(tran.start.date().isoformat(),'-')
    eventfile_datestring = datestr[0] + datestr[1] + datestr[2]
    
    event_list_dir=AFINO_EVENT_LIST_DIRECTORY #os.path.expanduser('~/afino_repository/event_lists/')
    if not os.path.exists(event_list_dir):
        os.makedirs(event_list_dir)
    pickfname=os.path.join(event_list_dir,'events_'+eventfile_datestring+'.pickle')

    #save the dictionary of events in a file for the given day - useful if we want to look over the events again later
    pickle.dump(events, open(pickfname, 'wb'))

    print(' ')
    print('---------------------------------------')
    print('Saved dictionary of events for ' + eventfile_datestring + ' to file: ' + pickfname)
    

    
        
def lightcurves_from_event_list(event_list_file, instrument='GBM'):
    '''This function reads an event list file and returns a lightcurve for each event in the file.
    The desired instrument may be selected with the 'instrument' keyword.'''

    print('Looking at list of lightcurves')
    data_dir = AFINO_DATA_DIRECTORY # os.path.expanduser('~/afino_repository/data/')
    event_list=pickle.load(open(event_list_file,'rb'))

    listoflightcurves=[]

    if len(event_list) == 0:
        print('No events in this event list. Returning')
        return []
        
    for event in event_list:
        start_date = event['start_time'].date().strftime('%Y%m%d')
        end_date = event['end_time'].date().strftime('%Y%m%d')
        goes_class = event['goes_class']

        start_date_fermi = event['start_time'].date().strftime('%y%m%d')
        allfiles=os.listdir(data_dir)

        #find the right file to open for the given instrument
        if instrument.upper() == 'GBM':
            #need to make edits here
            #file_to_open = [file for file in allfiles if file.startswith('glg_ctime') and (start_date_fermi in file)]
            file_to_open = '/Users/Inglis/afino_repository/detector_testing/glg_ctime_n3_121022_v00.pha'
            if not file_to_open:
                raise ValueError('Oops! No data files found in ' + data_dir +' matching the searched criteria')

            #find if there is a matching GBM trigger
            trig_date = datetime.datetime(event['start_time'].year,event['start_time'].month,event['start_time'].day)
            trig_start,trig_end = afino_gbm_utils.read_trigger_catalog(date = trig_date)
            #GBM trigger start is usually after GOES trigger
            #find events where GBM start trigger lies between GOES start and end times
            match_low = (trig_start > event['start_time'])
            match_high = (trig_start < event['end_time'])
            match = match_low * match_high
            gbm_trig_start = trig_start[match]
            gbm_trig_end = trig_end[match]
            stopping_time =  event['peak_time'] +  ((event['end_time'] - event['peak_time'])/2)
            print(gbm_trig_start)
            print(gbm_trig_end)
            print(stopping_time)
            
            if not gbm_trig_start.tolist():
                print('No matching events found in the GBM catalog. No lightcurve generated')
                lc2 = []
            else:
                print(file_to_open[0])
              #  lc2 = afino_gbm_utils.fermi_ctime_to_lightcurve_v2(os.path.join(data_dir,file_to_open[0]), gbm_trig_start[0], gbm_trig_end[-1], stopping_time)
                lc2 = afino_gbm_utils.fermi_ctime_to_lightcurve_v2(file_to_open, gbm_trig_start[0], gbm_trig_end[-1], stopping_time)
              #  lc = lightcurve.GBMSummaryLightCurve.create(file_to_open[0])
                lc_dframe = lc2.data.sort()
                lc2.data = lc_dframe
               # lc2 = lc.truncate(event['start_time'],event['end_time'])
                lc2.meta['GOES_CLASS'] = goes_class
        elif instrument.upper() == 'LYRA':
            file_to_open = [file for file in allfiles if file.startswith('lyra') and (start_date in file)]
            if not file_to_open:
                raise ValueError('Oops! No data files found in ' + data_dir +' matching the searched criteria')
            lc = lightcurve.LYRALightCurve.create(file_to_open[0])
            lc2 = lc.truncate(event['start_time'],event['end_time'])
            lc2.meta['GOES_CLASS'] = goes_class
        elif instrument.upper() == 'RHESSI':
            #rhessi a bit harder
            file_to_open=[file for file in allfiles if file.startswith('obs_summ')]
            if not file_to_open:
                raise ValueError('Oops! No data files found in ' + data_dir +' matching the searched criteria')
            lc = lightcurve.RHESSISummaryLightCurve.create(file_to_open[0])
            lc2 = lc.truncate(event['start_time'],event['end_time'])
            lc2.meta['GOES_CLASS'] = goes_class
        elif instrument.upper() == 'GOES':
            file_to_open=[file for file in allfiles if file.startswith('go') and (start_date in file)]
            if not file_to_open:
                raise ValueError('Oops! No data files found in ' + data_dir +' matching the searched criteria')
            lc = lightcurve.GOESLightCurve.create(file_to_open[0])
            lc2 = lc.truncate(event['start_time'],event['end_time'])
            lc2.meta['GOES_CLASS'] = goes_class
        else:
            raise ValueError('Instrument not recognised. Available options are GOES, GBM, RHESSI or LYRA')

        if lc2:
            listoflightcurves.append(lc2)

    print(listoflightcurves)
    return listoflightcurves
    

def prepare_lightcurves_for_analysis(listoflightcurves, wavelength, instrument='GBM'):
    '''Given a list of lightcurves, return a list of prepared TimeSeries that are ready to be input into the PyMC analysis.
    Returns a list of dictionaries - each dictionary contains a prepared signal for the chosen wavelength and instrument, plus some descriptor tags.

    Parameters
    ----------
    listoflightcurves - list
    A list of SunPy lightcurve objects, e.g. obtained for the flares of interest on a given date.

    wavelength - string
    A string descriptor of the wavelength of interest. Current allowed values are '4-15', '15-25', '25-50', '50-100' for GBM,
    and 'Al' or 'Zr' for LYRA.

    Keywords
    --------
    Instrument - string
    The instrument associated with the list of lightcurves. Current allowed values are 'GBM' and 'LYRA'.

    '''

    ref_time = parse_time('1979-01-01')
    #want to perform this operation on the list of lightcurves
    prepped_signals=[]
    for lc in listoflightcurves:
        #if the lightcurve contains no data then skip it
        if type(lc) == list:
            print('-------------------------------')
            print('Encountered a blank list in place of a lightcurve for ' + instrument + ' at ' + wavelength + '. Skipping preparation.')
            continue
        elif lc.data.empty:
            print('-------------------------------')
            print('Encountered a lightcurve with no data for ' + instrument + ' at ' + wavelength + '. Skipping preparation.')
            continue
        
        lc_dict={}
        times=[]
        for t in lc.data.index:
            times.append(parse_time(t.isoformat()))
        tt=[]
        for t in times:
            tt.append( (t - ref_time).total_seconds() )
        ttt = np.array(tt)

        allowed_gbm_wavelengths=['4-15','15-25','25-50','50-100']
        allowed_lyra_wavelengths=['Al','Zr']
        allowed_goes_wavelengths=['long']

        #exception handling             
        if (wavelength not in allowed_gbm_wavelengths) and (wavelength not in allowed_lyra_wavelengths) and (wavelength not in allowed_goes_wavelengths):
            raise ValueError('Specified wavelength for analysis is not valid. Returning. Acceptable values for wavelength parameter are: ' + allowed_gbm_wavelengths + ' ' + allowed_lyra_wavelengths + ' ' +
                             allowed_goes_wavelengths)
        elif instrument == 'GBM' and wavelength not in allowed_gbm_wavelengths:
            raise ValueError('Mismatch between specified instrument and analysis wavelength. Allowed GBM wavelengths are: ' + allowed_gbm_wavelengths)
        elif instrument == 'LYRA' and wavelength not in allowed_lyra_wavelengths:
            raise ValueError('Mismatch between specified instrument and analysis wavelength. Allowed LYRA wavelengths are: ' + allowed_lyra_wavelengths)
        elif instrument == 'GOES' and wavelength not in allowed_goes_wavelengths:
            raise ValueError('Mismatch between specified instrument and analysis wavelength. Allowed GOES wavelengths are: ' + allowed_goes_wavelengths)
        elif instrument == 'RHESSI':
            raise ValueError('RHESSI not currently supported. Returning.')

        #need to map the wavelength keys to the lc.data keys here. Use a dictionary.
        data_keys={}
        for wave in allowed_gbm_wavelengths:
            data_keys[wave] = wave+' keV'
        data_keys['Al'] = 'CHANNEL3'
        data_keys['Zr'] = 'CHANNEL4'
        data_keys['long'] = 'xrsb'

    
        #now create the TimeSeries object
        ts = TimeSeries(ttt,lc.data[data_keys[wavelength]])

        #prep the time series in different ways depending on the instrument. For GBM, have to filter out dropouts.
        #for LYRA, have to deal with LARs.
        #if instrument == 'GBM':
        #    ind2, therest = find_gbm_dropouts(lc)
        #    ts_prepped = treat_gbm_dropouts(ts,ind2,therest)
        #elif instrument == 'LYRA':
            #find and remove bad data regions (e.g. Large Angle Rotations)
        #    ts_prepped = ts
        #else:
        ts_prepped = ts

        lc_dict['timeseries'] = ts_prepped
        lc_dict['start_description'] = lc.time_range().start.strftime('%Y%m%d_%H%M%S')
        lc_dict['end_description'] = lc.time_range().end.strftime('%Y%m%d_%H%M%S')
        lc_dict['wavelength_description'] = wavelength
        lc_dict['instrument_description'] = instrument
        lc_dict['goes_class'] = lc.meta['GOES_CLASS']
        prepped_signals.append(lc_dict)
    
    return prepped_signals


    
        
    
def analyse_prepped_signals(prepped_signals,ppcheck=True):
    '''This routine takes in prepared lightcurve signals and performs the Bayesian MCMC analysis.'''
    for sig in prepped_signals:
        #first need to multiply by window function before input into the MCMC
        sig_apodized=prep_series(sig['timeseries'])
        #now perform model comparison
        description =  sig['start_description'] + '_' + sig['end_description'] + '_' + sig['instrument_description'] + '_' + sig['wavelength_description'] 
        root = '~/afino_repository/tempfiles/comparisons'
        model_comparison_dict=model_comparison(sig_apodized,root=root,description=description,ppcheck=ppcheck)

        #---------------------------------------------------------------------------------------------------------------------------
        #copy the files created by model comparison to somewhere more useful/consistent
        #this code block puts the model comparison summary files in ~/afino_repository/model_comparisons
        #MCMC analysis summary files are moved to ~/afino_repository/mcmc_analysis_results
        
        old_analysis_fnames=[os.path.join(os.path.expanduser(root),description,'analysis_summary.' + description + '_plaw_tmp.pickle'),
                             os.path.join(os.path.expanduser(root),description,'analysis_summary.' + description + '_gauss_tmp.pickle')]
        
        permanent_analysis_result_directory = PERMANENT_ANALYSIS_RESULT_DIRECTORY #os.path.expanduser('~/afino_repository/mcmc_analysis_results')
        final_analysis_fnames = [os.path.join(os.path.expanduser(permanent_analysis_result_directory),'analysis_summary.' + description + '_plaw.pickle'),
                                 os.path.join(os.path.expanduser(permanent_analysis_result_directory),'analysis_summary.' + description + '_gauss.pickle')]

        for i in range(0,len(final_analysis_fnames)):
            os.rename(old_analysis_fnames[i],final_analysis_fnames[i])

        old_model_comparison_fname = os.path.join(os.path.expanduser(root),description,description+'.pickle')
        permanent_model_comparison_summary_directory = MODEL_COMPARISON_DIRECTORY #os.path.expanduser('~/afino_repository/model_comparisons')
        final_model_comparison_fname = os.path.join(os.path.expanduser(permanent_model_comparison_summary_directory),'model_comparison_' + description + '.pickle')
        os.rename(old_model_comparison_fname,final_model_comparison_fname)
        
        #posterior predictive distribution files are moved to ~/afino_repository/posterior_predictive_results/
        old_posterior_predictive_fnames = [os.path.join(os.path.expanduser(root),description,'posterior_predictive.' + description + '_plaw_tmp.pickle'),
                                           os.path.join(os.path.expanduser(root),description,'posterior_predictive.' + description + '_gauss_tmp.pickle')]
        final_posterior_predictive_fnames = [os.path.join(POSTERIOR_PREDICTIVE_DIRECTORY,'posterior_predictive.' + description + '._plaw.pickle'),
                                             os.path.join(POSTERIOR_PREDICTIVE_DIRECTORY,'posterior_predictive.' + description + '._gauss.pickle')]

        for i in range(0,len(final_posterior_predictive_fnames)):
            os.rename(old_posterior_predictive_fnames[i], final_posterior_predictive_fnames[i])
        #------------------------------------------------------------------------------------------------------------------------------

        if model_comparison_dict['BIC_ratio'] > 0.0:
            best_model='splwc_AddNormalBump2'
            p_values = model_comparison_dict['p_values_gauss']
            final_posterior_predictive_fname = final_posterior_predictive_fnames[1]
        else:
            best_model='single_power_law_with_constant'
            p_values = model_comparison_dict['p_values_plaw']
            final_posterior_predictive_fname = final_posterior_predictive_fnames[0]
            
     #   if ppcheck:    
      #      pcheck_results = main_analysis(sig_apodized,ppcheck=True,plots=False,root=root,description=description,
       #                   nsample=1000,model=best_model)
            #more to follow here
        #p_values = pcheck_results[4]
            
        create_analysis_record(sig,final_analysis_fnames,final_model_comparison_fname, p_values = p_values)
        create_analysis_summary_plot(sig, final_analysis_fnames, final_model_comparison_fname, final_posterior_predictive_fname = final_posterior_predictive_fname, p_values = p_values)


def analyse_prepped_signals2(prepped_signals,ppcheck=True):
    '''This routine takes in prepared lightcurve signals and performs the Bayesian MCMC analysis.'''
    for sig in prepped_signals:
        #first need to multiply by window function before input into the MCMC
        sig_apodized=prep_series(sig['timeseries'])
        #now perform model comparison
        description =  sig['start_description'] + '_' + sig['end_description'] + '_' + sig['instrument_description'] + '_' + sig['wavelength_description'] 
        root = '~/afino_repository/tempfiles/comparisons'
        analysis_summary=model_comparison(sig_apodized,root=root,description=description)


        if analysis_summary['dBIC'] > 0.0:
            best_model='splwc_AddNormalBump2'
        else:
            best_model='single_power_law_with_constant'
            
        create_analysis_record2(sig,analysis_summary)
     #   create_analysis_summary_plot(sig, final_analysis_fnames, final_model_comparison_fname, final_posterior_predictive_fname = final_posterior_predictive_fname, p_values = p_values)
        create_analysis_summary_plot3(sig, analysis_summary)
        
def create_analysis_record(sig, analysis_filenames, model_comparison_filename, p_values = None):
    '''
    This routine creates an analysis record saved to a python pickle file. This record
    contains a summary of the Bayesian MCMC analysis of an event.
    '''

    analysis_record_directory = ANALYSIS_RECORD_DIRECTORY #os.path.expanduser('~/afino_repository/analysis_records')
    description =  sig['start_description'] + '_' + sig['end_description'] + '_' + sig['instrument_description'] + '_' + sig['wavelength_description'] 
    analysis_record={}
    res = pickle.load(open(model_comparison_filename,'rb'))

    analysis_record['event'] = sig['start_description'] + '_' + sig['end_description']
    analysis_record['instrument'] = sig['instrument_description']
    analysis_record['wavelength'] = sig['wavelength_description']
    analysis_record['goes_class'] = sig['goes_class']
    analysis_record['bic'] = res['BIC_ratio']
    analysis_record['model1'] = 'single_power_law_plus_constant'
    analysis_record['model2'] = 'spwlc_normalbump'
    analysis_record['model_comparison_file'] = model_comparison_filename
    analysis_record['analysis_filename'] = analysis_filenames

    #also want to save the MAP values here, not the means
    if (res['BIC_ratio'] < 0.0):
        params = pickle.load(open(analysis_filenames[0],'rb'))
        #analysis_record['fit_params'] = [params['stats']['power_law_index']['mean']]
        analysis_record['fit_params'] = [params['map_power_law_index']]
    else:
        params = pickle.load(open(analysis_filenames[1],'rb'))
       # analysis_record['fit_params'] = [params['stats']['power_law_index']['mean'],params['stats']['gaussian_position']['mean'],params['stats']['gaussian_width']['mean']]
        analysis_record['fit_params'] = [params['map_power_law_index'],params['map_gaussian_position'],params['map_gaussian_width']]

    analysis_record['period'] = []
    if len(analysis_record['fit_params']) > 1:
        analysis_record['period'] = 1 / (np.exp(analysis_record['fit_params'][1]))

     #if available, add pcheck results to record
    if p_values:
        analysis_record['p_value_tsse'] = round(p_values['p_value_tsse'],5)
        analysis_record['p_value_tr'] = round(p_values['p_value_tr'],5)

    #add some additional flags to the analysis record
    # S = short time series (< 200 data points)
    # B = bad fit to preferred model
    analysis_record['flags'] = []
    if len(sig['timeseries'].data) < 200:
        analysis_record['flags'] +='S'
    if (p_values['p_value_tsse'] < 0.01) or (p_values['p_value_tsse'] > 0.99):
        analysis_record['flags'] +='B'

    

        
    
    pickle.dump(analysis_record,open(os.path.join(analysis_record_directory,'analysis_record_'+description+'.pkl'),'wb'))


def create_analysis_record2(sig, analysis_summary):
    '''
    This routine creates an analysis record saved to a python pickle file. This record
    contains a summary of the Bayesian MCMC analysis of an event.
    '''

    analysis_record_directory = ANALYSIS_RECORD_DIRECTORY #os.path.expanduser('~/afino_repository/analysis_records')
    description =  sig['start_description'] + '_' + sig['end_description'] + '_' + sig['instrument_description'] + '_' + sig['wavelength_description'] 
    analysis_record={}
   # res = pickle.load(open(model_comparison_filename,'rb'))

    analysis_record['event'] = sig['start_description'] + '_' + sig['end_description']
    analysis_record['instrument'] = sig['instrument_description']
    analysis_record['wavelength'] = sig['wavelength_description']
    analysis_record['goes_class'] = sig['goes_class']
    analysis_record['bic'] = analysis_summary['dBIC']
    analysis_record['dBIC_0v2'] = analysis_summary['dBIC_0v2']
    analysis_record['dBIC_2v1'] = analysis_summary['dBIC_2v1']
    analysis_record['model0'] = 'single_power_law_plus_constant'
    analysis_record['model1'] = 'spwlc_AddNormalBump2'
    analysis_record['model2'] = 'broken_power_law_with_constant'
    #analysis_record['model_comparison_file'] = model_comparison_filename
    analysis_record['analysis_filename'] = analysis_summary['analysis_filename']
    
    #also want to save the best fit params for both models:
    analysis_record['m0_fit_params'] = analysis_summary['m0']['params']
    analysis_record['m1_fit_params'] = analysis_summary['m1']['params']
    analysis_record['m2_fit_params'] = analysis_summary['m2']['params']

    #also save the goodness of fit stats for both models.
    analysis_record['rchi2_m0'] = analysis_summary['m0']['rchi2']
    analysis_record['probability_m0'] = analysis_summary['m0']['probability']
    analysis_record['rchi2_m1'] = analysis_summary['m1']['rchi2']
    analysis_record['probability_m1'] = analysis_summary['m1']['probability']
    analysis_record['rchi2_m2'] = analysis_summary['m2']['rchi2']
    analysis_record['probability_m2'] = analysis_summary['m2']['probability']
    

    #also want to save the MAP values here, not the means
  #  if (res['BIC_ratio'] < 0.0):
   #     params = pickle.load(open(analysis_filenames[0],'rb'))
        #analysis_record['fit_params'] = [params['stats']['power_law_index']['mean']]
    #    analysis_record['fit_params'] = [params['map_power_law_index']]
    #else:
     #   params = pickle.load(open(analysis_filenames[1],'rb'))
       # analysis_record['fit_params'] = [params['stats']['power_law_index']['mean'],params['stats']['gaussian_position']['mean'],params['stats']['gaussian_width']['mean']]
      #  analysis_record['fit_params'] = [params['map_power_law_index'],params['map_gaussian_position'],params['map_gaussian_width']]

   # analysis_record['period'] = []
    if analysis_summary['dBIC'] > 0.0:
        analysis_record['period'] = 1 / (np.exp(analysis_record['m1_fit_params'][4]))
        
        
     #if available, add pcheck results to record
#    if p_values:
 #       analysis_record['p_value_tsse'] = round(p_values['p_value_tsse'],5)
  #      analysis_record['p_value_tr'] = round(p_values['p_value_tr'],5)

    #add some additional flags to the analysis record
    # S = short time series (< 200 data points)
    # B0, B1 = bad fit to model 0, 1
    analysis_record['flags'] = []
    if len(sig['timeseries'].data) < 200:
        analysis_record['flags'] +='S'

    if analysis_record['probability_m0'] < 0.01:
        analysis_record['flags'].append('B0')
    if analysis_record['probability_m1'] < 0.01:
        analysis_record['flags'].append('B1')
    if analysis_record['probability_m2'] < 0.01:
        analysis_record['flags'].append('B2')

    

        
    
    pickle.dump(analysis_record,open(os.path.join(analysis_record_directory,'analysis_record_'+description+'.pkl'),'wb'))
    return
        
def create_master_record(gbm_only=False,goes_only=True):
    '''
    This function creates a master record from all of the individual analysis records, and saves
    the master record in a CSV file.
    '''

    
    analysis_record_directory = ANALYSIS_RECORD_DIRECTORY #os.path.expanduser('~/afino_repository/analysis_records')
    csvfile = open(os.path.join(analysis_record_directory,'afino_master_record.csv'), 'wb')
    #allfiles = os.listdir(analysis_record_directory)
    records = glob.glob(analysis_record_directory+'/*.pkl')
    fieldnames = ['Date','GOES_class','Start_time','End_time','Instrument','Wavelength','dBIC_0v1','dBIC_0v2','dBIC_2v1','Detection','rchi2_m0','probability_m0',
                  'rchi2_m1','probability_m1','rchi2_m2','probability_m2','period','width','Flags']
    writer = csv.DictWriter(csvfile,fieldnames=fieldnames)
    writer.writeheader()
    c = 0
    for record in records:
        print(c)
        c +=1
        result = pickle.load(open(os.path.join(analysis_record_directory,record),'rb'))
        if gbm_only:
            if result['instrument'] != 'GBM':
                continue
        elif goes_only:
            if result['instrument'] != 'GOES':
                continue
                
        result_dict = {}
        date_info = result['event'].split('_')
        result_dict['Date'] = date_info[0]
        result_dict['GOES_class'] = result['goes_class']
        result_dict['Start_time'] = date_info[1]
        result_dict['End_time'] = date_info[3]
        result_dict['Instrument'] = result['instrument']
        if result['wavelength'] == 'long':
            result_dict['Wavelength'] = '1-8A'
        else:
            result_dict['Wavelength'] = result['wavelength']
        
        result_dict['dBIC_0v1'] = str(result['bic'])[0:8]
        result_dict['dBIC_0v2'] = str(result['dBIC_0v2'])[0:8]
        result_dict['dBIC_2v1'] = str(result['dBIC_2v1'])[0:8]
        
        if (result['bic'] > 10.0) and (result['dBIC_2v1'] > 10.0):
            result_dict['Detection'] = 'Yes'
        else:
            result_dict['Detection'] = 'No'

        #add the goodness-of-fit values
        if 'rchi2_m0' in result:
            result_dict['rchi2_m0'] = result['rchi2_m0']
            result_dict['probability_m0'] = result['probability_m0']
            result_dict['rchi2_m1'] = result['rchi2_m1']
            result_dict['probability_m1'] = result['probability_m1']
            result_dict['rchi2_m2'] = result['rchi2_m2']
            result_dict['probability_m2'] = result['probability_m2']
            
        if 'flags' in result:
            result_dict['Flags'] = result['flags']
        else:
            result_dict['Flags'] = []


        #if there was a detection, store the period
        if 'period' in result:
            result_dict['period'] = round(result['period'],1)
            result_dict['width'] = result['m1_fit_params'][5]

        if gbm_only:
            #look through the GBM bad event list and add a flag if there is a match
            bad_event_list = csv.DictReader(open('/Users/Inglis/physics/proba2_analysis/gbm_bad_mask3.csv','rU'))
            bad = False
            for ev in bad_event_list:
                if (ev['Date'] == date_info[0] and ev['Start time'] == date_info[1]):
                    bad = True
            if bad:
                result_dict['Flags'].append('D')
        
        writer.writerow(result_dict)
    return



def create_analysis_summary_plot(sig, analysis_filenames, model_comparison_filename, final_posterior_predictive_fname = None, p_values = None):

    #directory paths
    model_comparison_directory=MODEL_COMPARISON_DIRECTORY #os.path.expanduser('~/afino_repository/model_comparisons')
    permanent_analysis_result_directory = PERMANENT_ANALYSIS_RESULT_DIRECTORY #os.path.expanduser('~/afino_repository/mcmc_analysis_results')
    analysis_summary_plots_directory = ANALYSIS_SUMMARY_PLOTS_DIRECTORY #os.path.expanduser('~/afino_repository/analysis_summary_plots')
    
    
    import seaborn as sns
    sns.set_style('ticks',{'xtick.direction':'in','ytick.direction':'in'})
    sns.set_context('paper')

    ts = sig['timeseries']
    basetime = parse_time(ts.SampleTimes.basetime)

    model_comparison = pickle.load(open(model_comparison_filename,'rb'))
    BIC = model_comparison['BIC_ratio']
    
    plt.figure(1,figsize=(12,10))
    plt.subplots_adjust(bottom=0.067,top=0.967,left=0.1,right=0.95)
    
    plt.subplot(2,2,1)
    plt.tick_params(labelsize=12)
    plt.plot(ts.SampleTimes.time,ts.data)

    btime2=basetime.isoformat()
    btimestr=btime2[0:19]
    plt.xlabel('Start time ' + btimestr)
    plt.ylabel('Intensity (arb.)')
    plt.title(sig['instrument_description'] + ' ' + sig['wavelength_description'],fontsize=16)
    #plt.xlim([100,500])

    plt.subplot(2,2,2)
    plt.tick_params(labelsize=12)

    summ_file = os.path.join(permanent_analysis_result_directory, analysis_filenames[0])
    summ=pickle.load(open(summ_file,'rb'))

    #this is wrong - need to plot the MAP best-fit params instead of the mean. Necessary params are in 'summ'
    plt.loglog(summ['frequencies'],summ['power'],label='data')
    #plt.loglog(summ['frequencies'],summ['stats']['fourier_power_spectrum']['mean'],label='best fit')
    plt.loglog(summ['frequencies'],summ['map_fourier_power_spectrum'],label='best_fit')
    plt.legend(fontsize=14)
    plt.xlabel('frequency (Hz)',fontsize=12)
    plt.ylabel('Fourier power',fontsize=12)
    plt.title('Power Spectral Density (PSD)',fontsize=12)

    plt.figtext(0.61,0.65,r'$\alpha=%4.2f$' % summ['stats']['power_law_index']['mean'],fontsize=12)
    plt.figtext(0.61,0.63,r'$A=%4.3f$' % summ['stats']['power_law_norm']['mean'],fontsize=12)
    

    plt.subplot(2,2,3)
    if final_posterior_predictive_fname:
        post = pickle.load(open(final_posterior_predictive_fname,'rb'))
        tsse_distribution = post[1]['vaughan_2010_T_SSE']
        observed_tsse = post[0]['vaughan_2010_T_SSE']
        plt.hist(tsse_distribution,bins=30)
        plt.axvline(observed_tsse)
        plt.figtext(0.38,0.4,'p = ' + str(round(p_values['p_value_tsse'],5) ))
        plt.title('Test statistic $T_{sse}$ distribution')
    else:
        plt.plot([])
        plt.xlim(0,1)
        ply.ylim(0,1)
        plt.title('Test statistic $T_{sse}$ distribution')
        plt.figtext(0.35,0.8,'p = Null')
    plt.figtext(0.38,0.37,'$\Delta$BIC = ' +str(round(BIC,3)))
    
    plt.subplot(2,2,4)
    plt.tick_params(labelsize=12)

    summ_file = os.path.join(permanent_analysis_result_directory, analysis_filenames[1])
    summ=pickle.load(open(summ_file,'rb'))
    
    plt.loglog(summ['frequencies'],summ['power'],label='data')
    #plt.loglog(summ['frequencies'],summ['stats']['fourier_power_spectrum']['mean'],label='best fit')
    plt.loglog(summ['frequencies'],summ['map_fourier_power_spectrum'],label='best_fit')
    plt.legend(fontsize=14)
    plt.xlabel('frequency (Hz)',fontsize=12)
    plt.ylabel('Fourier power',fontsize=12)
    #plt.title('Power Spectral Density (PSD)',fontsize=12)

    period=1/np.exp(summ['stats']['gaussian_position']['mean'])
    plt.axvline(np.exp(summ['stats']['gaussian_position']['mean']),color='red',linestyle='--')
    plt.figtext(0.61,0.17,r'$\alpha=%4.2f$' % summ['stats']['power_law_index']['mean'],fontsize=12)
    plt.figtext(0.61,0.15,r'$A=%4.3f$' % np.exp(summ['stats']['gaussian_amplitude']['mean']),fontsize=12)
    plt.figtext(0.61,0.13,r'$f_0=%4.3f$' % np.exp(summ['stats']['gaussian_position']['mean']) + ' Hz = '
                + r'$%4.2f$' % period +' s',fontsize=12)
    plt.figtext(0.61,0.11,r'$\sigma=%4.3f$' % summ['stats']['gaussian_width']['mean'],fontsize=12)

    savefilename = 'summary_plot_' + sig['start_description'] + '_' + sig['end_description'] + '_' + sig['instrument_description'] + '_' + sig['wavelength_description'] + '.pdf'
    plt.savefig(os.path.join(analysis_summary_plots_directory, savefilename ))
    plt.close()
    


def create_analysis_summary_plot2(sig, analysis_summary):

    #directory paths
    model_comparison_directory=MODEL_COMPARISON_DIRECTORY #os.path.expanduser('~/afino_repository/model_comparisons')
    permanent_analysis_result_directory = PERMANENT_ANALYSIS_RESULT_DIRECTORY #os.path.expanduser('~/afino_repository/mcmc_analysis_results')
    analysis_summary_plots_directory = ANALYSIS_SUMMARY_PLOTS_DIRECTORY #os.path.expanduser('~/afino_repository/analysis_summary_plots')
    
    
    import seaborn as sns
    sns.set_style('ticks',{'xtick.direction':'in','ytick.direction':'in'})
    sns.set_context('paper')

    ts = sig['timeseries']
    basetime = parse_time(ts.SampleTimes.basetime)

   # model_comparison = pickle.load(open(model_comparison_filename,'rb'))
    BIC = analysis_summary['dBIC'] #model_comparison['BIC_ratio']
    
    plt.figure(1,figsize=(18,5))
    plt.subplots_adjust(bottom=0.1,top=0.9,left=0.05,right=0.95)
    
    plt.subplot(1,3,1)
    plt.tick_params(labelsize=12)
    plt.plot(ts.SampleTimes.time,ts.data)

    btime2=basetime.isoformat()
    btimestr=btime2[0:19]
    plt.xlabel('Start time ' + btimestr)
    plt.ylabel('Intensity (arb.)')
    if sig['wavelength_description'] == 'long':
        plt.title('GOES 1-8$\AA$',fontsize=16)
    else:
        plt.title(sig['instrument_description'] + ' ' + sig['wavelength_description'],fontsize=16)
    #plt.xlim([100,500])

    plt.subplot(1,3,2)
    plt.tick_params(labelsize=12)

   # summ_file = os.path.join(permanent_analysis_result_directory, analysis_filenames[0])
    #summ=pickle.load(open(summ_file,'rb'))

    #this is wrong - need to plot the MAP best-fit params instead of the mean. Necessary params are in 'summ'
    plt.loglog(analysis_summary['m0']['frequencies'],analysis_summary['m0']['power'],label='data')
    #plt.loglog(summ['frequencies'],summ['stats']['fourier_power_spectrum']['mean'],label='best fit')
    plt.loglog(analysis_summary['m0']['frequencies'],analysis_summary['m0']['best_fit_power_spectrum'],label='best fit')
    plt.legend(fontsize=14)
    plt.xlabel('frequency (Hz)',fontsize=12)
    plt.ylabel('Fourier power',fontsize=12)
    plt.title('Power Spectral Density (PSD) - Model 1',fontsize=12)

    plt.figtext(0.40,0.25,r'$\alpha=%4.2f$' % analysis_summary['m0']['params'][1],fontsize=12)
    plt.figtext(0.58,0.70,r'$\chi^2=%4.2f$' % analysis_summary['m0']['rchi2']) 
   # plt.figtext(0.61,0.63,r'$A=%4.3f$' % analysis_summary['m0']['params'][0],fontsize=12)
    

    plt.figtext(0.23,0.77,'$\Delta$BIC = ' +str(round(BIC,3)))
    
    plt.subplot(1,3,3)
    plt.tick_params(labelsize=12)

  #  summ_file = os.path.join(permanent_analysis_result_directory, analysis_filenames[1])
   # summ=pickle.load(open(summ_file,'rb'))
    
    plt.loglog(analysis_summary['m1']['frequencies'],analysis_summary['m1']['power'],label='data')
   
    plt.loglog(analysis_summary['m1']['frequencies'],analysis_summary['m1']['best_fit_power_spectrum'],label='best fit')
    plt.legend(fontsize=14)
    plt.xlabel('frequency (Hz)',fontsize=12)
    plt.ylabel('Fourier power',fontsize=12)
    plt.title('Power Spectral Density (PSD) - Model 2',fontsize=12)

    period=1/np.exp(analysis_summary['m1']['params'][4])
    plt.axvline(np.exp(analysis_summary['m1']['params'][4]),color='red',linestyle='--')
    plt.figtext(0.72,0.25,r'$\alpha=%4.2f$' % analysis_summary['m1']['params'][1],fontsize=12)
  #  plt.figtext(0.61,0.15,r'$A=%4.3f$' % np.exp(analysis_summary['stats']['gaussian_amplitude']['mean']),fontsize=12)
    plt.figtext(0.72,0.21,r'$f_0=%4.3f$' % np.exp(analysis_summary['m1']['params'][4]) + ' Hz = '
                + r'$%4.2f$' % period +' s',fontsize=12)
    plt.figtext(0.72,0.17,r'$\sigma=%4.3f$' % analysis_summary['m1']['params'][5],fontsize=12)
    plt.figtext(0.9,0.70,r'$\chi^2=%4.2f$' % analysis_summary['m1']['rchi2']) 

    savefilename = 'summary_plot_' + sig['start_description'] + '_' + sig['end_description'] + '_' + sig['instrument_description'] + '_' + sig['wavelength_description'] + '.pdf'
    plt.savefig(os.path.join(analysis_summary_plots_directory, savefilename ))
    plt.close()

    return

        

def create_analysis_summary_plot3(sig, analysis_summary):

    #directory paths
    model_comparison_directory=MODEL_COMPARISON_DIRECTORY #os.path.expanduser('~/afino_repository/model_comparisons')
    permanent_analysis_result_directory = PERMANENT_ANALYSIS_RESULT_DIRECTORY #os.path.expanduser('~/afino_repository/mcmc_analysis_results')
    analysis_summary_plots_directory = ANALYSIS_SUMMARY_PLOTS_DIRECTORY #os.path.expanduser('~/afino_repository/analysis_summary_plots')
    
    
    import seaborn as sns
    sns.set_style('ticks',{'xtick.direction':'in','ytick.direction':'in'})
    sns.set_context('paper')

    ts = sig['timeseries']
    basetime = parse_time(ts.SampleTimes.basetime)

   # model_comparison = pickle.load(open(model_comparison_filename,'rb'))
    BIC = analysis_summary['dBIC'] #model_comparison['BIC_ratio']
    
    plt.figure(1,figsize=(24,5))
    plt.subplots_adjust(bottom=0.1,top=0.9,left=0.05,right=0.95)
    
    plt.subplot(1,4,1)
    plt.tick_params(labelsize=12)
    plt.plot(ts.SampleTimes.time,ts.data)

    btime2=basetime.isoformat()
    btimestr=btime2[0:19]
    plt.xlabel('Start time ' + btimestr,fontsize=14)
    plt.ylabel('Intensity (arb.)')
    if sig['wavelength_description'] == 'long':
        plt.title('GOES 1-8$\AA$',fontsize=16)
    else:
        plt.title(sig['instrument_description'] + ' ' + sig['wavelength_description'],fontsize=16)
    #plt.xlim([100,500])

    plt.subplot(1,4,2)
    plt.tick_params(labelsize=12)

   # summ_file = os.path.join(permanent_analysis_result_directory, analysis_filenames[0])
    #summ=pickle.load(open(summ_file,'rb'))

    #this is wrong - need to plot the MAP best-fit params instead of the mean. Necessary params are in 'summ'
    plt.loglog(analysis_summary['m0']['frequencies'],analysis_summary['m0']['power'],label='data')
    #plt.loglog(summ['frequencies'],summ['stats']['fourier_power_spectrum']['mean'],label='best fit')
    plt.loglog(analysis_summary['m0']['frequencies'],analysis_summary['m0']['best_fit_power_spectrum'],label='best fit')
    plt.legend(fontsize=14)
    plt.xlabel('frequency (Hz)',fontsize=12)
    plt.ylabel('Fourier power',fontsize=12)
    plt.title('Power Spectral Density (PSD) - Model 0',fontsize=12)

    plt.figtext(0.31,0.25,r'$\alpha=%4.2f$' % analysis_summary['m0']['params'][1],fontsize=14)
    plt.figtext(0.43,0.50,r'$\chi^2=%4.2f$' % analysis_summary['m0']['rchi2'],fontsize=14) 
   # plt.figtext(0.61,0.63,r'$A=%4.3f$' % analysis_summary['m0']['params'][0],fontsize=12)
    

    plt.figtext(0.165,0.77,'$\Delta$BIC$_{M0-M1}$ = ' +str(round(BIC,1)),fontsize=14)
    plt.figtext(0.165,0.72,'$\Delta$BIC$_{M0-M2}$ = ' + str(round(analysis_summary['dBIC_0v2'],1)),fontsize=14)
  

    
    plt.subplot(1,4,3)
    plt.tick_params(labelsize=12)


    #work out the 1sig and 2sig upper levels for model m1
    m1_plaw = rnspectralmodels3.power_law_with_constant(analysis_summary['m1']['params'][0:3],analysis_summary['m1']['frequencies'])
    m1_quantile_1sig = (-np.log(0.16) * m1_plaw)
    m1_quantile_2sig = (-np.log(0.025) * m1_plaw)
    m1_quantile_2sig_lower = (-np.log(0.975) * m1_plaw)

    
  #  summ_file = os.path.join(permanent_analysis_result_directory, analysis_filenames[1])
   # summ=pickle.load(open(summ_file,'rb'))
    
    plt.loglog(analysis_summary['m1']['frequencies'],analysis_summary['m1']['power'],label='data')
   
    plt.loglog(analysis_summary['m1']['frequencies'],analysis_summary['m1']['best_fit_power_spectrum'],label='best fit')
    plt.loglog(analysis_summary['m1']['frequencies'],m1_quantile_2sig,color='orangered',linestyle='dashed',label='2$\sigma$')
    plt.loglog(analysis_summary['m1']['frequencies'],m1_quantile_2sig_lower,color='orangered',linestyle='dashed')
    plt.legend(fontsize=14)
    plt.xlabel('frequency (Hz)',fontsize=12)
    plt.ylabel('Fourier power',fontsize=12)
    plt.title('Power Spectral Density (PSD) - Model 1',fontsize=12)

    period=1/np.exp(analysis_summary['m1']['params'][4])
    plt.axvline(np.exp(analysis_summary['m1']['params'][4]),color='red',linestyle='--')
    plt.figtext(0.54,0.25,r'$\alpha=%4.2f$' % analysis_summary['m1']['params'][1],fontsize=14)
  #  plt.figtext(0.61,0.15,r'$A=%4.3f$' % np.exp(analysis_summary['stats']['gaussian_amplitude']['mean']),fontsize=12)
    plt.figtext(0.54,0.21,r'$f_0=%4.3f$' % np.exp(analysis_summary['m1']['params'][4]) + ' Hz = '
                + r'$%4.2f$' % period +' s',fontsize=14)
    plt.figtext(0.54,0.17,r'$\sigma=%4.3f$' % analysis_summary['m1']['params'][5],fontsize=14)
    plt.figtext(0.67,0.5,r'$\chi^2=%4.2f$' % analysis_summary['m1']['rchi2'],fontsize=14) 


    plt.subplot(1,4,4)
    plt.tick_params(labelsize=12)

   # summ_file = os.path.join(permanent_analysis_result_directory, analysis_filenames[0])
    #summ=pickle.load(open(summ_file,'rb'))

    #this is wrong - need to plot the MAP best-fit params instead of the mean. Necessary params are in 'summ'
    plt.loglog(analysis_summary['m2']['frequencies'],analysis_summary['m2']['power'],label='data')
    #plt.loglog(summ['frequencies'],summ['stats']['fourier_power_spectrum']['mean'],label='best fit')
    plt.loglog(analysis_summary['m2']['frequencies'],analysis_summary['m2']['best_fit_power_spectrum'],label='best fit')
    plt.legend(fontsize=14)
    plt.xlabel('frequency (Hz)',fontsize=12)
    plt.ylabel('Fourier power',fontsize=12)
    plt.title('Power Spectral Density (PSD) - Model 2',fontsize=12)

    plt.figtext(0.78,0.25,r'$\alpha=%4.2f$' % analysis_summary['m2']['params'][1],fontsize=14)
    plt.figtext(0.78,0.21,r'$\alpha_2=%4.2f$' % analysis_summary['m2']['params'][3],fontsize=14)
    plt.figtext(0.90,0.5,r'$\chi^2=%4.2f$' % analysis_summary['m2']['rchi2'],fontsize=14)
    
   # plt.figtext(0.61,0.63,r'$A=%4.3f$' % analysis_summary['m0']['params'][0],fontsize=12)
    




    
    savefilename = 'summary_plot_' + sig['start_description'] + '_' + sig['end_description'] + '_' + sig['instrument_description'] + '_' + sig['wavelength_description'] + '.pdf'
    plt.savefig(os.path.join(analysis_summary_plots_directory, savefilename ))
    plt.close()

    return

        

        
        
        
        
        
        
    

    

