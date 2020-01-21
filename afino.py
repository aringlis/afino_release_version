
# AFINO - Automated Flare INference of Oscillations
#-----------------------------------------------------------------
#
# AUTHOR : Andrew Inglis
#
# PURPOSE : Perform automated Bayesian MCMC of solar flare lightcurves from various instruments.
# The analysis is based on the procedure presented in Inglis et al., ApJ, 2015.
# The goal is to identify events of interest that may contain an oscillatory component in their
# temporal signal. 
#
# PROCEDURE : Following Inglis et al. (2015, 2016), the AFINO code consists of the following steps:
#
#             1) Search GOES event list for a given date for solar flare events
#             2) Download data for any events found for GOES, GBM and LYRA instruments
#             3) Prepare lightcurves for each instrument and wavelength during the event timeframe
#             4) Pass the prepared lightcurves into the Bayesian MCMC analysis routine
#             5) Create an analysis record summarizing the results of the event analysis
#             6) Create a summary plot for each successfully completed event analysis
#             7) Store the model comparison and fit results for future reference
#             8) Update the master record for every event analysed
#
# DATA PRODUCTS : A master record in CSV format is kept showing the results of every event analysed
#                 Additionally, each fit and model comparison is permanently stored, along with a
#                 analysis record summary file and plot for each event.
#
# UPDATED : 28th October 2015, ARI
#

#this is the main AFINO executable
import os
import afino_utils

def main(date):
    #some presets
    event_list_dir = os.path.expanduser('~/afino_repository/event_lists/')
    instruments = ['GOES']#,'GBM']#,'LYRA']
    gbm_wavelengths = ['15-25']
    lyra_wavelengths = ['Al','Zr']
    goes_wavelengths = ['long']
    
    #search the GOES catalogue for solar flares on the given date
    afino_utils.search_for_events_on_date(date)
    event_file = os.path.join(event_list_dir,'events_'+ date.strftime('%Y%m%d') + '.pickle')



    for instr in instruments:
        #generate a list of lightcurves from the constructed daily event file for chosen instrument
        try:
            listoflightcurves=afino_utils.lightcurves_from_event_list(event_file,instrument=instr)
        except ValueError as e:
            print(e)
            continue
        if len(listoflightcurves) == 0:
            print('Empty lightcurve list encountered. Skipping')
            continue
    
        #then prepare the signal for analysis by the Bayesian MCMC method and run the analysis
        if instr == 'GBM':
            print('-------------------')
            print('Analysing data from ' + instr)
            print('-------------------')
            for wave in gbm_wavelengths:
                print('---------------')
                print('Analysing signal for ' + instr + ' at ' + wave)
                print('---------------')
                prepped_signals = afino_utils.prepare_lightcurves_for_analysis(listoflightcurves,wave,instrument=instr)
                #run the Bayesian MCMC analysis
                try:
                    afino_utils.analyse_prepped_signals2(prepped_signals)
                except RuntimeError:
                    print('RuntimeError occured while trying to perform Bayesian MCMC model comparison for ' + instr + ' at ' + wave)
        elif instr == 'LYRA':
            print('-------------------')
            print('Analysing data from ' + instr)
            print('-------------------')
            for wave in lyra_wavelengths:
                print('---------------')
                print('Analysing signal for ' + instr + ' at ' + wave)
                print('---------------')
                prepped_signals = afino_utils.prepare_lightcurves_for_analysis(listoflightcurves,wave,instrument=instr)
                #run the Bayesian MCMC analysis
                try:
                    afino_utils.analyse_prepped_signals2(prepped_signals)
                except RuntimeError:
                    print('RuntimeError occured while trying to perform Bayesian MCMC model comparison for ' + instr + ' at ' + wave)

        elif instr == 'GOES':
            print('-------------------')
            print('Analysing data from ' + instr)
            print('-------------------')
            for wave in goes_wavelengths:
                print('---------------')
                print('Analysing signal for ' + instr + ' at ' + wave)
                print('---------------')
                prepped_signals = afino_utils.prepare_lightcurves_for_analysis(listoflightcurves,wave,instrument=instr)
                #run the Bayesian MCMC analysis
                try:
                    afino_utils.analyse_prepped_signals2(prepped_signals)
                except RuntimeError:
                    print('RuntimeError occured while trying to perform Bayesian MCMC model comparison for ' + instr + ' at ' + wave)

   # afino_utils.create_master_record()


    
