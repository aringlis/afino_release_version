import matplotlib.pyplot as plt
import prettyplotlib as ppl
import pickle
import numpy as np
from sunpy.time import parse_time
import os


def compare_noise_levels(summ_file,pymc_file):
    summ=pickle.load(open(summ_file,'rb'))
    pymc=pickle.load(open(pymc_file,'rb'))


    #plot some simulated data on top
    #pick some random simulations
    l=len(pymc['predictive'][0])
    num=np.random.randint(0,high=l,size=2)
    
    plt.loglog(summ['frequencies'],pymc['predictive'][0][num[0]],label='simulated data')
    plt.loglog(summ['frequencies'],pymc['predictive'][0][num[1]],label='simulated data')
    #plt.loglog(summ['frequencies'],pymc['predictive'][0][num[2]])

    #plot the original Fourier spectrum
    plt.loglog(summ['frequencies'],summ['power'],label='data')
    #best fit signal
    plt.loglog(summ['frequencies'],summ['stats']['fourier_power_spectrum']['mean'],label='best fit')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('Fourier power')
    
    plt.legend()
    plt.show()
    
def plot_fourier_spectrum(file,post_file):

    summ=pickle.load(open(file,'rb'))
    post=pickle.load(open(post_file,'rb'))
    #original signal
    plt.loglog(summ['frequencies'],summ['power'],label='data')
    #best fit signal
    plt.loglog(summ['frequencies'],summ['stats']['fourier_power_spectrum']['mean'],label='best fit power law,' + r'$\alpha=%4.2f$' % summ['stats']['power_law_index']['mean'])

    #T_R is MAX(R) where R = 2 Iobs / Sj
    #Sj is the model spectrum
    #Iobs is the observed spectrum
    #So want to find the value of Iobs that would result in a certain T_R
    #Iobs = Sj x R /2 where R is the T_R value corresponding to a certain level, e.g. 99%
    #First step is to find R from the distribution
    ordered=np.sort(post[1]['vaughan_2010_T_R'])
    #put the T_Rs in ascending order
    l=len(ordered)
    #find the R corresponding to the 99% level
    r99=ordered[l*0.99]
    #now calculate the 99% Fourier spectral density line based on that
    l99 = summ['stats']['fourier_power_spectrum']['mean'] * r99 / 2    

    #total hack follows
    #l99= summ['stats']['fourier_power_spectrum']['mean'] * 22. / 2
    #l95= summ['stats']['fourier_power_spectrum']['mean'] * 18. / 2
    #l66= summ['stats']['fourier_power_spectrum']['mean'] * 14.2 / 2

    plt.loglog(summ['frequencies'], l99, label='99% confidence level')
    #plt.loglog(summ['frequencies'], l95, label='95% confidence level')
    #plt.loglog(summ['frequencies'], l66, label='66% confidence level')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('Fourier power')
    plt.legend()
    plt.show()

def plot_distribution(file):

    r=pickle.load(open(file,'rb'))
    a=r[0]
    dtr=r[1]['vaughan_2010_T_R']
    dtsse=r[1]['vaughan_2010_T_SSE']
    
    plt.figure(1)
    h=plt.hist(dtr,bins=20)
    plt.axvline(a['vaughan_2010_T_R'],color='red',linewidth=2)
    plt.xlabel('Statistic value')
    plt.ylabel('N')
    plt.title('Statistic: Vaughan 2010 T_R')
    plt.show()

    plt.figure(2)
    h2=plt.hist(dtsse,bins=20)
    plt.axvline(a['vaughan_2010_T_SSE'],color='red',linewidth=2)
    plt.xlabel('Statistic value')
    plt.ylabel('N')
    plt.title('Statistic: Vaughan 2010 T_SSE')
    plt.show()

def plot_original_series(file):

    summ=pickle.load(open(file,'rb'))
    plt.loglog(summ['frequencies'],summ['power'])

def plot_all(i,ts,summ_file,post_file,sub=None,file_desc='default',savedir=None):

    plt.figure(1,figsize=(12,12))
    plt.subplot(2,2,1)
    plt.plot(ts.SampleTimes.time,ts.data / np.hanning(len(ts.data)))
    btime=parse_time(np.double(ts.SampleTimes.basetime))
    btime2=btime.isoformat()
    btimestr=btime2[0:19]
    plt.xlabel('Start time ' + btimestr)
    plt.ylabel('Flux (normalised)')
    #plot the right title
    if file_desc == 'l3':
        if sub:
            plt.title('LYRA Al filter '+'('+str(sub)+')')
        else:
            plt.title('LYRA Al filter')
    elif file_desc == 'l4':
        if sub:
            plt.title('LYRA Zr filter '+'('+str(sub)+')')
        else:
            plt.title('LYRA Zr filter')
    elif file_desc == '612':
        plt.title('RHESSI 6-12 keV')
    elif file_desc == '1225':
        plt.title('RHESSI 12-25 keV')
    elif file_desc == '2550':
        plt.title('RHESSI 25-50 keV')
    elif file_desc == '50100':
        plt.title('RHESSI 50-100 keV')
    else:
        plt.title('Time Series')
    


    summ=pickle.load(open(summ_file,'rb'))
    r=pickle.load(open(post_file,'rb'))
    #original signal
    plt.subplot(2,2,2)
    plt.loglog(summ['frequencies'],summ['power'],label='data')
    #best fit signal
    plt.loglog(summ['frequencies'],summ['stats']['fourier_power_spectrum']['mean'],label='best fit power law,' + r'$\alpha=%4.2f$' % summ['stats']['power_law_index']['mean'])

    #T_R is MAX(R) where R = 2 Iobs / Sj
    #Sj is the model spectrum
    #Iobs is the observed spectrum
    #So want to find the value of Iobs that would result in a certain T_R
    #Iobs = Sj x R /2 where R is the T_R value corresponding to a certain level, e.g. 99%
    #First step is to find R from the distribution
    ordered=np.sort(r[1]['vaughan_2010_T_R'])
    #put the T_Rs in ascending order
    l=len(ordered)
    #find the R corresponding to the 99% level
    r99=ordered[l*0.99]
    #now calculate the 99% Fourier spectral density line based on that
    l99 = summ['stats']['fourier_power_spectrum']['mean'] * r99 / 2    

    plt.loglog(summ['frequencies'], l99, label='99% confidence level')
    plt.axvline(x=0.0625,color='red')
    plt.xlabel('frequency (Hz)')
    plt.ylabel('Fourier power')
    plt.legend(framealpha=0)


	
    
    a=r[0]
    dtr=r[1]['vaughan_2010_T_R']
    dtsse=r[1]['vaughan_2010_T_SSE']

    #work out the p-values for the distributions
    p_tr=len(np.where(dtr > r[0]['vaughan_2010_T_R'])[0]) / np.float(len(dtr))
    p_tsse=len(np.where(dtsse > r[0]['vaughan_2010_T_SSE'])[0]) / np.float(len(dtsse))

    
    plt.subplot(2,2,3)
    h=plt.hist(dtr,bins=20)
    plt.axvline(a['vaughan_2010_T_R'],color='red',linewidth=2)
    plt.xlabel('Statistic value')
    plt.ylabel('N')
    plt.figtext(0.35,0.37,'p = '+str(p_tr),fontsize=16)
    plt.title('Statistic: Vaughan 2010 T_R')

    plt.subplot(2,2,4)
    h2=plt.hist(dtsse,bins=20)
    plt.axvline(a['vaughan_2010_T_SSE'],color='red',linewidth=2)
    plt.xlabel('Statistic value')
    plt.ylabel('N')
    plt.figtext(0.78,0.37,'p = '+str(p_tsse),fontsize=16)
    plt.title('Statistic: Vaughan 2010 T_SSE')
    if sub:
        plt.savefig(os.path.join(savedir,'summary_fig_event_'+file_desc+'_'+str(i)+'_subinterval_'+str(sub)+'.pdf'))
    else:
        plt.savefig(os.path.join(savedir,'summary_fig_event_'+file_desc+'_'+str(i)+'.pdf'))
    plt.close()
    

def plot_paper_figure(i,ts,summ_file,post_file,sub=None,file_desc='default',savedir=None,fourpanel=False):

    import seaborn as sns
    sns.set_style('ticks',{'xtick.direction':'in','ytick.direction':'in'})
    sns.set_context('paper')

    if fourpanel:
        plt.figure(1,figsize=(24,5))
    else:
        plt.figure(1,figsize=(12,4))  
        
    plt.subplots_adjust(bottom=0.15,left=0.085,right=0.95)

    if fourpanel:
        plt.subplot(1,4,1)
    else:
        plt.subplot(1,2,1)

    plt.tick_params(labelsize=14)
    plt.plot(ts.SampleTimes.time,ts.data)#  / np.hanning(len(ts.data)))
    btime=parse_time(np.double(ts.SampleTimes.basetime))
    btime2=btime.isoformat()
    btimestr=btime2[0:19]
    plt.xlabel('Time (s) from ' + btimestr,fontsize=14)
    plt.ylabel('Flux (normalised)',fontsize=14)
    
    #plot the right title

    if file_desc == 'l3':
        if sub:
            plt.title('LYRA Al filter '+'('+str(sub)+')',fontsize=18)
            if sub == 1:
                plt.figtext(0.095,0.8,'e)',fontsize=18)
        else:
            plt.title('LYRA Al filter',fontsize=18)
            plt.figtext(0.095,0.8,'d)',fontsize=18)
    elif file_desc == 'l4':
        if sub:
            plt.title('LYRA Zr filter '+'('+str(sub)+')',fontsize=18)
            if sub == 1:
                plt.figtext(0.095,0.8,'g)',fontsize=18)
        else:
            plt.title('LYRA Zr filter',fontsize=18)
            plt.figtext(0.095,0.8,'f)',fontsize=18)
    elif file_desc == '612':
        plt.title('RHESSI 6-12 keV',fontsize=18)
    elif file_desc == '1225':
        plt.title('RHESSI 12-25 keV',fontsize=18)
    elif file_desc == '2550':
        plt.title('RHESSI 25-50 keV',fontsize=18)
    elif file_desc == '50100':
        plt.title('RHESSI 50-100 keV',fontsize=18)
    elif file_desc == 'g415' or file_desc == 'g415_n2':
        plt.title('GBM 4-15 keV',fontsize=18)
        plt.figtext(0.095,0.8,'a)',fontsize=18)
    elif file_desc == 'g1227' or file_desc == 'g1227_n2':
        plt.title('GBM 12-27 keV',fontsize=18)
        plt.figtext(0.095,0.8,'a)',fontsize=18)
    elif file_desc == 'g2750' or file_desc == 'g2750_n2':
        plt.title('GBM 27-50 keV',fontsize=18)
        plt.figtext(0.095,0.8,'b)',fontsize=18)
    elif file_desc == 'g50100' or file_desc == 'g50100_n2':
        plt.title('GBM 50-100 keV',fontsize=18)
        plt.figtext(0.095,0.8,'c)',fontsize=18)
    else:
        plt.title('Time Series',fontsize=18)
    


    summ=pickle.load(open(summ_file,'rb'))
    r=pickle.load(open(post_file,'rb'))
    #original signal
    if fourpanel:
        plt.subplot(1,4,2)
    else:
        plt.subplot(1,2,2)

    plt.tick_params(labelsize=14)
    plt.loglog(summ['frequencies'],summ['power'],label='data')
    #best fit signal
    plt.loglog(summ['frequencies'],summ['stats']['fourier_power_spectrum']['mean'],label='best fit, ' + r'$\alpha=%4.2f$' % summ['stats']['power_law_index']['mean'])
    plt.xlim([5e-4,5e-1])
    plt.title('Power spectral density (PSD)',fontsize=18)
    #T_R is MAX(R) where R = 2 Iobs / Sj
    #Sj is the model spectrum
    #Iobs is the observed spectrum
    #So want to find the value of Iobs that would result in a certain T_R
    #Iobs = Sj x R /2 where R is the T_R value corresponding to a certain level, e.g. 99%
    #First step is to find R from the distribution
    ordered=np.sort(r[1]['vaughan_2010_T_R'])
    #put the T_Rs in ascending order
    l=len(ordered)
    #find the R corresponding to the 99% level
    r99=ordered[l*0.99]
    #now calculate the 99% Fourier spectral density line based on that
    l99 = summ['stats']['fourier_power_spectrum']['mean'] * r99 / 2    

    plt.loglog(summ['frequencies'], l99, label='99% interval')
    #plt.axvline(x=0.0625,color='red')
    plt.xlabel('frequency (Hz)',fontsize=14)
    plt.ylabel('Fourier power',fontsize=14)
    plt.legend(framealpha=0,fontsize=16)

    
    a=r[0]
    dtr=r[1]['vaughan_2010_T_R']
    dtsse=r[1]['vaughan_2010_T_SSE']

    #work out the p-values for the distributions
    p_tr=len(np.where(dtr > r[0]['vaughan_2010_T_R'])[0]) / np.float(len(dtr))
    p_tsse=len(np.where(dtsse > r[0]['vaughan_2010_T_SSE'])[0]) / np.float(len(dtsse))

    if fourpanel:
        plt.subplot(1,4,3)
        plt.tick_params(labelsize=14)
        h=plt.hist(dtr,bins=20)
        plt.axvline(a['vaughan_2010_T_R'],color='red',linewidth=2)
        plt.xlabel('Statistic value',fontsize=14)
        plt.ylabel('N',fontsize=14)
        plt.figtext(0.65,0.77,'p = '+str(p_tr),fontsize=16)
        plt.title('$T_R$ distribution',fontsize=18)

        plt.subplot(1,4,4)
        plt.tick_params(labelsize=14)
        h2=plt.hist(dtsse,bins=20)
        plt.axvline(a['vaughan_2010_T_SSE'],color='red',linewidth=2)
        plt.xlabel('Statistic value',fontsize=14)
        plt.ylabel('N',fontsize=14)
        plt.figtext(0.89,0.77,'p = '+str(p_tsse),fontsize=16)
        plt.title('$T_{SSE}$ distribution',fontsize=18)

        if file_desc in ['l3','l4'] and sub is not None:
            plt.savefig(os.path.join(savedir,'paper_fig_event_'+file_desc+'_'+str(i)+'_subinterval_'+str(sub)+'4panel.pdf'))
        else:
            plt.savefig(os.path.join(savedir,'paper_fig_event_'+file_desc+'_'+str(i)+'4panel.pdf'))
            
    else:
        if file_desc in ['l3','l4'] and sub is not None:
            plt.savefig(os.path.join(savedir,'paper_fig_event_'+file_desc+'_'+str(i)+'_subinterval_'+str(sub)+'.pdf'))
        else:
            plt.savefig(os.path.join(savedir,'paper_fig_event_'+file_desc+'_'+str(i)+'.pdf'))
        
    plt.close()
    
        
