
import os
import datetime
import numpy as np

from afino.afino_series import AfinoSeries, prep_series
from afino.afino_model_comparison import model_comparison, model_comparison_mms
from afino import afino_spectral_models
import matplotlib.pyplot as plt


def analyse_series(times, flux, description=None, low_frequency_cutoff=None, high_frequency_cutoff = None, savedir=None, overwrite_gauss_bounds = None,
                       overwrite_extra_gauss_bounds = None, use_json = True, model_ids = [0,1,2]):
    """
    Analyse a single, generic timeseries using the AFINO model comparison code.
    
    Parameters
    ----------

    times : ndarray
        an array of times
    flux : ndarray
        an array of data points
    description : string, optional
        a string descriptor of the analysis run, that is incorporated into output filenames
    low_frequency_cutoff : float, optional
        specifies a frequency above which the input Fourier spectrum is not analysed
    savedir : string, optional
        specifies a directory for output save files
    use_json : bool
        If set to True, saves analysis output in JSON format. If False, pickle format is used.
        Default is True.
    
    """

    ts = AfinoSeries(times,flux)
    #first need to apply a window function 
    sig_apodized=prep_series(ts)
    #now perform model comparison
    if not description:
        description = datetime.datetime.now().strftime('%Y%m%d_%H%M%S') 
   
    analysis_summary=model_comparison(sig_apodized,description=description, low_frequency_cutoff=low_frequency_cutoff, high_frequency_cutoff = high_frequency_cutoff,
                                    overwrite_gauss_bounds = overwrite_gauss_bounds, overwrite_extra_gauss_bounds = overwrite_extra_gauss_bounds, use_json = use_json,
                                    model_ids = model_ids)
 
    
   # create_generic_summary_plot(ts, analysis_summary, description, low_frequency_cutoff=low_frequency_cutoff, savedir=savedir)

    return analysis_summary

def analyse_series_twobump(times, flux, description=None, low_frequency_cutoff=None, savedir=None, overwrite_gauss_bounds = None, overwrite_extra_gauss_bounds = None):
    """Analyse a single, generic timeseries using the AFINO model comparison code., with extra models"""

    ts = AfinoSeries(times,flux)
    #first need to multiply by window function before input into the analysis
    sig_apodized=prep_series(ts)
    #now perform model comparison
    if not description:
        description = 'afino_series_' + datetime.datetime.now().strftime('%Y%m%d_%H%M%S') 
    root = '~/afino_repository/tempfiles/comparisons'
   
    analysis_summary=model_comparison_mms(sig_apodized,root=root,description=description, generic=True, low_frequency_cutoff=low_frequency_cutoff,
                                    overwrite_gauss_bounds = overwrite_gauss_bounds , overwrite_extra_gauss_bounds = overwrite_extra_gauss_bounds)

  
    
    create_generic_summary_plot_twobump(ts, analysis_summary, description, low_frequency_cutoff=low_frequency_cutoff, savedir=savedir)
    
    
    
def create_generic_summary_plot(ts, analysis_summary, description, low_frequency_cutoff=None, savedir=None):
    """Create a summary plot showing the result of an AFINO analysis run."""
    
    import seaborn as sns
    sns.set_style('ticks',{'xtick.direction':'in','ytick.direction':'in'})
    sns.set_context('paper')

    #ensure the directory to save plots exists, create it if not.
    os.makedirs(os.path.expanduser('~/afino_repository/plots/'),exist_ok=True)

    npanels = len(analysis_summary) + 1

    plt.figure(1,figsize=(6*npanels,5))
    plt.subplots_adjust(bottom=0.1,top=0.9,left=0.05,right=0.95)

    plt.subplot(1,npanels,1)
    plt.tick_params(labelsize=12)
    plt.plot(ts.SampleTimes.time,ts.data)

    plt.xlabel('Time (s)',fontsize=14)
    plt.ylabel('Intensity (arb.)')
    plt.title('Input signal',fontsize=16)
    #plt.xlim([100,500])

    
   # for i in range(0,npanels-1):
    for i, key in enumerate(analysis_summary.keys()):
        #key = 'm' + str(i)
        plt.subplot(1,npanels,i+2)
        plt.tick_params(labelsize=12)

        plt.loglog(analysis_summary[key]['frequencies'],analysis_summary[key]['power'],label='data')
        plt.loglog(analysis_summary[key]['frequencies'],analysis_summary[key]['best_fit_power_spectrum'],label='best fit')
        plt.legend(fontsize=14)
        plt.xlabel('frequency (Hz)',fontsize=12)
        plt.ylabel('Fourier power',fontsize=12)
        plt.title('Power Spectral Density (PSD) - Model ' + str(analysis_summary[key]['ID']),fontsize=12)

        plt.xlim([1e-5,1e2])
        plt.text(2e-5,1e-2,r'$\alpha=%4.2f$' % analysis_summary[key]['params'][1],fontsize=14)
        plt.text(2e-3,1.0,r'$\chi^2=%4.2f$' % analysis_summary[key]['rchi2'],fontsize=14) 

     #   plt.figtext(0.165,0.77,'$\Delta$BIC$_{M0-M1}$ = ' +str(round(BIC,1)),fontsize=14)
      #  plt.figtext(0.165,0.72,'$\Delta$BIC$_{M0-M2}$ = ' + str(round(analysis_summary['dBIC_0v2'],1)),fontsize=14)

        if low_frequency_cutoff:
            plt.axvline(low_frequency_cutoff)

        if (analysis_summary[key]['model'] == 'pow_const_gauss') or (analysis_summary[key]['model'] == 'pow_const_2gauss'):
            #work out the 1sig and 2sig upper levels for model m1
            m1_plaw = afino_spectral_models.pow_const(analysis_summary[key]['params'][0:3],analysis_summary[key]['frequencies'])
            m1_quantile_1sig = (-np.log(0.16) * m1_plaw)
            m1_quantile_2sig = (-np.log(0.025) * m1_plaw)
            m1_quantile_2sig_lower = (-np.log(0.975) * m1_plaw)
            plt.loglog(analysis_summary[key]['frequencies'],m1_quantile_2sig,color='orangered',
                           linestyle='dashed',label='2$\sigma$')
            plt.loglog(analysis_summary[key]['frequencies'],m1_quantile_2sig_lower,color='orangered',
                           linestyle='dashed')

            # display the best-fit period location
            period=1/np.exp(analysis_summary[key]['params'][4])
            plt.axvline(np.exp(analysis_summary[key]['params'][4]),color='red',linestyle='--')
            plt.text(2e-5,1e-3,r'$f_0=%4.3f$' % np.exp(analysis_summary[key]['params'][4]) + ' Hz = '
                + r'$%4.2f$' % period +' s',fontsize=14)

            if analysis_summary[key]['model'] == 'pow_const_2gauss':
                period2=1/np.exp(analysis_summary[key]['params'][7])
                plt.axvline(np.exp(analysis_summary[key]['params'][7]),color='red',linestyle='--')
                plt.text(2e-5,1e-4,r'$f_1=%4.3f$' % np.exp(analysis_summary[key]['params'][7]) + ' Hz = '
                + r'$%4.2f$' % period2 +' s',fontsize=14)
            
        if analysis_summary[key]['model'] == 'bpow_const':
            plt.text(2e-5,1e-3,r'$\alpha_2=%4.2f$' % analysis_summary[key]['params'][3],fontsize=14)

 
            
            

    

    savefilename = 'afino_summary_plot_' + description + '.pdf'
    if savedir:
        plt.savefig(os.path.join(savedir,savefilename))
    else:
        plt.savefig(os.path.join(os.path.expanduser('~/afino_repository/plots/'), savefilename ))
    plt.close()

    return



def create_generic_summary_plot_twobump(ts, analysis_summary, description, low_frequency_cutoff=None, savedir=None):
    """Create a summary plot showing the result of an AFINO analysis run."""
    
    import seaborn as sns
    sns.set_style('ticks',{'xtick.direction':'in','ytick.direction':'in'})
    sns.set_context('paper')

    #ensure the directory to save plots exists, create it if not.
    os.makedirs(os.path.expanduser('~/afino_repository/plots/'),exist_ok=True)
    
   # model_comparison = pickle.load(open(model_comparison_filename,'rb'))
    BIC = analysis_summary['dBIC'] #model_comparison['BIC_ratio']

    #get relative BIC values for easy comparison
    bic_list = [analysis_summary['m0']['BIC'], analysis_summary['m1']['BIC'], analysis_summary['m2']['BIC'], analysis_summary['m3']['BIC'], analysis_summary['m4']['BIC']]
    loc = np.argmin(bic_list)
    maxloc = np.argmax(bic_list)
    bic_min = bic_list[loc]
    bic_max = bic_list[maxloc]

    
    
    plt.figure(1,figsize=(36,5))
    plt.subplots_adjust(bottom=0.1,top=0.9,left=0.05,right=0.95)
    
    plt.subplot(1,6,1)
    plt.tick_params(labelsize=12)
    plt.plot(ts.SampleTimes.time,ts.data)

    plt.xlabel('Time (s)',fontsize=14)
    plt.ylabel('Intensity (arb.)')
    plt.title('Generic Series',fontsize=16)
    #plt.xlim([100,500])

    plt.subplot(1,6,2)
    plt.tick_params(labelsize=12)

   # summ_file = os.path.join(permanent_analysis_result_directory, analysis_filenames[0])
    #summ=pickle.load(open(summ_file,'rb'))

    plt.loglog(analysis_summary['m0']['frequencies'],analysis_summary['m0']['power'],label='data')
    plt.loglog(analysis_summary['m0']['frequencies'],analysis_summary['m0']['best_fit_power_spectrum'],label='best fit')
    plt.legend(fontsize=14)
    plt.xlabel('frequency (Hz)',fontsize=12)
    plt.ylabel('Fourier power',fontsize=12)
    plt.title('Power Spectral Density (PSD) - Model 0',fontsize=12)

    plt.figtext(0.22,0.25,r'$\alpha=%4.2f$' % analysis_summary['m0']['params'][1],fontsize=14)
    plt.figtext(0.30,0.65,r'$\chi^2=%4.2f$' % analysis_summary['m0']['rchi2'],fontsize=14) 

    plt.figtext(0.14,0.77,'$BIC_0$ = ' + str(round(bic_list[0] - bic_max,1)),fontsize=14,alpha=0.8)
    plt.figtext(0.14,0.72,'$BIC_1$ = ' + str(round(bic_list[1] - bic_max,1)),fontsize=14,alpha=0.8)
    plt.figtext(0.14,0.67,'$BIC_2$ = ' + str(round(bic_list[2] - bic_max,1)),fontsize=14,alpha=0.8)
    plt.figtext(0.14,0.62,'$BIC_3$ = ' + str(round(bic_list[3] - bic_max,1)),fontsize=14,alpha=0.8)
    plt.figtext(0.14,0.57,'$BIC_4$ = ' + str(round(bic_list[4] - bic_max,1)),fontsize=14,alpha=0.8)

    
    if low_frequency_cutoff:
        plt.axvline(low_frequency_cutoff)
    
    plt.subplot(1,6,3)
    plt.tick_params(labelsize=12)

    #work out the 1sig and 2sig upper levels for model m1
    m1_plaw = afino_spectral_models.power_law_with_constant(analysis_summary['m1']['params'][0:3],analysis_summary['m1']['frequencies'])
    m1_quantile_1sig = (-np.log(0.16) * m1_plaw)
    m1_quantile_2sig = (-np.log(0.025) * m1_plaw)
    m1_quantile_2sig_lower = (-np.log(0.975) * m1_plaw)


    
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
    plt.figtext(0.38,0.25,r'$\alpha=%4.2f$' % analysis_summary['m1']['params'][1],fontsize=14)
  #  plt.figtext(0.61,0.15,r'$A=%4.3f$' % np.exp(analysis_summary['stats']['gaussian_amplitude']['mean']),fontsize=12)
    plt.figtext(0.38,0.21,r'$f_0=%4.3f$' % np.exp(analysis_summary['m1']['params'][4]) + ' Hz = '
                + r'$%4.2f$' % period +' s',fontsize=14)
    plt.figtext(0.38,0.17,r'$\sigma=%4.3f$' % analysis_summary['m1']['params'][5],fontsize=14)
    plt.figtext(0.45,0.65,r'$\chi^2=%4.2f$' % analysis_summary['m1']['rchi2'],fontsize=14) 

    if low_frequency_cutoff:
        plt.axvline(low_frequency_cutoff)

    plt.subplot(1,6,4)
    plt.tick_params(labelsize=12)

   # summ_file = os.path.join(permanent_analysis_result_directory, analysis_filenames[0])
    #summ=pickle.load(open(summ_file,'rb'))

    plt.loglog(analysis_summary['m2']['frequencies'],analysis_summary['m2']['power'],label='data')
    plt.loglog(analysis_summary['m2']['frequencies'],analysis_summary['m2']['best_fit_power_spectrum'],label='best fit')
    plt.legend(fontsize=14)
    plt.xlabel('frequency (Hz)',fontsize=12)
    plt.ylabel('Fourier power',fontsize=12)
    plt.title('Power Spectral Density (PSD) - Model 2',fontsize=12)

    plt.figtext(0.52,0.25,r'$\alpha=%4.2f$' % analysis_summary['m2']['params'][1],fontsize=14)
    plt.figtext(0.52,0.21,r'$\alpha_2=%4.2f$' % analysis_summary['m2']['params'][3],fontsize=14)
    plt.figtext(0.60,0.65,r'$\chi^2=%4.2f$' % analysis_summary['m2']['rchi2'],fontsize=14)

    if low_frequency_cutoff:
        plt.axvline(low_frequency_cutoff)
   # plt.figtext(0.61,0.63,r'$A=%4.3f$' % analysis_summary['m0']['params'][0],fontsize=12)


    plt.subplot(1,6,5)
    plt.tick_params(labelsize=12)

   # summ_file = os.path.join(permanent_analysis_result_directory, analysis_filenames[0])
    #summ=pickle.load(open(summ_file,'rb'))

    plt.loglog(analysis_summary['m3']['frequencies'],analysis_summary['m3']['power'],label='data')
    plt.loglog(analysis_summary['m3']['frequencies'],analysis_summary['m3']['best_fit_power_spectrum'],label='best fit')
    plt.legend(fontsize=14)
    plt.xlabel('frequency (Hz)',fontsize=12)
    plt.ylabel('Fourier power',fontsize=12)
    plt.title('Power Spectral Density (PSD) - Model 3',fontsize=12)

    period=1/np.exp(analysis_summary['m3']['params'][4])
    plt.axvline(np.exp(analysis_summary['m3']['params'][4]),color='red',linestyle='--')
    plt.figtext(0.68,0.35,r'$\alpha=%4.2f$' % analysis_summary['m1']['params'][1],fontsize=14)
  #  plt.figtext(0.61,0.15,r'$A=%4.3f$' % np.exp(analysis_summary['stats']['gaussian_amplitude']['mean']),fontsize=12)
    plt.figtext(0.68,0.30,r'$f_0=%4.3f$' % np.exp(analysis_summary['m3']['params'][4]) + ' Hz = '
                + r'$%4.2f$' % period +' s',fontsize=14)
    plt.figtext(0.68,0.25,r'$\sigma=%4.3f$' % analysis_summary['m3']['params'][5],fontsize=14)
    plt.figtext(0.76,0.65,r'$\chi^2=%4.2f$' % analysis_summary['m3']['rchi2'],fontsize=14) 

    period2 = 1/np.exp(analysis_summary['m3']['params'][7])
    plt.axvline(np.exp(analysis_summary['m3']['params'][7]),color='red',linestyle='--')
    plt.figtext(0.68,0.20,r'$f_{0b}=%4.3f$' % np.exp(analysis_summary['m3']['params'][7]) + ' Hz = '
                + r'$%4.2f$' % period2 +' s',fontsize=14)
    plt.figtext(0.68,0.15,r'$\sigma_b=%4.3f$' % analysis_summary['m3']['params'][8],fontsize=14)

    
    if low_frequency_cutoff:
        plt.axvline(low_frequency_cutoff)


    plt.subplot(1,6,6)

    plt.loglog(analysis_summary['m4']['frequencies'],analysis_summary['m4']['power'],label='data')
   
    plt.loglog(analysis_summary['m4']['frequencies'],analysis_summary['m4']['best_fit_power_spectrum'],label='best fit')
    plt.loglog(analysis_summary['m4']['frequencies'],m1_quantile_2sig,color='orangered',linestyle='dashed',label='2$\sigma$')
    plt.loglog(analysis_summary['m4']['frequencies'],m1_quantile_2sig_lower,color='orangered',linestyle='dashed')
    plt.legend(fontsize=14)
    plt.xlabel('frequency (Hz)',fontsize=12)
    plt.ylabel('Fourier power',fontsize=12)
    plt.title('Power Spectral Density (PSD) - Model 4',fontsize=12)

    period3=1/np.exp(analysis_summary['m4']['params'][4])
    plt.axvline(np.exp(analysis_summary['m4']['params'][4]),color='red',linestyle='--')
    plt.figtext(0.85,0.25,r'$\alpha=%4.2f$' % analysis_summary['m4']['params'][1],fontsize=14)
  #  plt.figtext(0.61,0.15,r'$A=%4.3f$' % np.exp(analysis_summary['stats']['gaussian_amplitude']['mean']),fontsize=12)
    plt.figtext(0.85,0.21,r'$f_0=%4.3f$' % np.exp(analysis_summary['m4']['params'][4]) + ' Hz = '
                + r'$%4.2f$' % period3 +' s',fontsize=14)
    plt.figtext(0.85,0.17,r'$\sigma=%4.3f$' % analysis_summary['m4']['params'][5],fontsize=14)
    plt.figtext(0.91,0.65,r'$\chi^2=%4.2f$' % analysis_summary['m4']['rchi2'],fontsize=14) 

    if low_frequency_cutoff:
        plt.axvline(low_frequency_cutoff)
    

    savefilename = 'summary_plot_' + description + '.pdf'
    if savedir:
        plt.savefig(os.path.join(savedir,savefilename))
    else:
        plt.savefig(os.path.join(os.path.expanduser('~/afino_repository/plots/'), savefilename ))
    plt.close()

    return
