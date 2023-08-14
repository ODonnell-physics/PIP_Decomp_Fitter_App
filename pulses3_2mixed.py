# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 10:57:39 2022

@author: Brian
"""
import os 
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import fftconvolve as fftconvolve

def model(t, a1, mu1, sig1, tau1, a2, mu2, sig2, tau2, a3, mu3, sig3, tau3, w):
    try:
        measured_signal2
    except:
        measured_signal2 = [0, 0]

    pulse1 = a1 * np.exp(-.5 * ((t - mu1) / sig1) ** 2)
    pulse2 = a2 * np.exp(-.5 * ((t - mu2) / sig2) ** 2)

    decay_fast = np.exp(-(t) / tau1)
    decay_slow = np.exp(-(t) / tau2)

    sumdecay = (w * decay_fast) + (1-w)*decay_slow
    norm_sumdecay = sumdecay / max(sumdecay)

    flux1 = fftconvolve(pulse1, norm_sumdecay)
    flux2 = fftconvolve(pulse2, norm_sumdecay)

    pulse3 = a3 * np.exp(-.5 * ((t - mu3) / sig3) ** 2)
    decay3 = np.exp(-(t) / tau3)

    flux3 = fftconvolve(pulse3, decay3)

    test_signal = (flux1 + flux2 + flux3) + np.min(np.abs(measured_signal2))

    t_l = len(t)

    test_signal = (flux1 + flux2 + flux3)[0:t_l] + np.min(np.abs(measured_signal2))

    return test_signal

def p3_2m(m_s,m_s_e,freqs,p_times, f_i,f_f,t_i,t_f,ig,ub,lb, plotter=True, f_unit="GHz", b_unit="SFU", t_unit="s", log_scale=False, rolling_ig=False):
    
    def model(t,a1,mu1,sig1,tau1,a2,mu2,sig2,tau2,a3,mu3,sig3,tau3,w):
    
        pulse1=a1*np.exp(-.5*((t-mu1)/sig1)**2)
        pulse2=a2*np.exp(-.5*((t-mu2)/sig2)**2)
            
        decay_fast=np.exp(-(t)/tau1)
        decay_slow=np.exp(-(t)/tau2)
            
        sumdecay = (w * decay_fast) + (1-w)*decay_slow
        norm_sumdecay=sumdecay/max(sumdecay)
        
        flux1=fftconvolve(pulse1, norm_sumdecay)
        flux2=fftconvolve(pulse2, norm_sumdecay)
    
        
        pulse3=a3*np.exp(-.5*((t -mu3)/sig3)**2)
        decay3=np.exp(-(t)/tau3)
        
        flux3=fftconvolve(pulse3, decay3)

        t_l= len(t)
        
        test_signal=(flux1+flux2+flux3)[0:t_l]+np.min(np.abs(measured_signal2))
    
        return test_signal
    
    plt.rcParams['font.size'] = 24
    

    
    

    

    f_w= f_f - f_i
    s_l= t_f - t_i
    
    
    
    result=np.zeros([13,f_w])
    result_error=np.zeros([13,f_w])
    fl=np.around(freqs,4)
    
    
    xirs=np.zeros([f_w])

    s_l_=s_l
    if s_l % 2 != 0: s_l_ = s_l + 1
    pulse1_spec = np.zeros([f_w, s_l_])
    pulse2_spec = np.zeros([f_w, s_l_])
    pulse3_spec = np.zeros([f_w, s_l_])

    flux1_spec = np.zeros([f_w, s_l_])
    flux2_spec = np.zeros([f_w, s_l_])
    flux3_spec = np.zeros([f_w, s_l_])

    oddflag= False

    if plotter:
        if not os.path.exists('./3_pulses_2_mixed_fit'):
            os.makedirs('./3_pulses_2_mixed_fit')

        pref = './3_pulses_2_mixed_fit/3_2'



    for f in range(f_i, f_f, 1):

        measured_signal2 = m_s[t_i:t_f, f]
        measured_signal2e = m_s_e[t_i:t_f, f]


        time= p_times[t_i:t_f]-p_times[t_i]

        if s_l % 2 != 0:
            s_l = s_l + 1
            oddflag=True

        if s_l != len(measured_signal2):
            measured_signal2 = np.append(measured_signal2, m_s[t_f, f])
            measured_signal2e = np.append(measured_signal2e, m_s[t_f, f])
            time= np.append(time,time[s_l-2])


        l = 1.0



        try:
            fit = curve_fit(model, time, measured_signal2[0:s_l], sigma=measured_signal2e[0:s_l], p0=ig, absolute_sigma=True, bounds=(lb, ub))


            ans,cov = fit
        
            fit_a1,fit_mu1,fit_sig1,fit_tau1,fit_a2,fit_mu2,fit_sig2,fit_tau2,fit_a3,fit_mu3,fit_sig3,fit_tau3,fit_w = ans
            fit_sa1,fit_smu1,fit_ssig1,fit_stau1,fit_sa2,fit_smu2,fit_ssig2,fit_stau2,fit_sa3,fit_smu3,fit_ssig3,fit_stau3,fit_sw  = np.sqrt(np.diag(cov))

            if rolling_ig: ig = ans
        
            results=np.array([fit_a1,fit_mu1,fit_sig1,fit_tau1,fit_a2,fit_mu2,fit_sig2,fit_tau2,fit_a3,fit_mu3,fit_sig3,fit_tau3,fit_w ])
            results_error= np.array([fit_sa1,fit_smu1,fit_ssig1,fit_stau1,fit_sa2,fit_smu2,fit_ssig2,fit_stau2,fit_sa3,fit_smu3,fit_ssig3,fit_stau3,fit_sw])
        
            
                
            pulse1=fit_a1*np.exp(-.5*((time-fit_mu1)/fit_sig1)**2)
            pulse2=fit_a2*np.exp(-.5*((time-fit_mu2)/fit_sig2)**2)
                
            decay_fast=np.exp(-(time)/fit_tau1)
            decay_slow=np.exp(-(time)/fit_tau2)
                
            sumdecay= fit_w*decay_fast+decay_slow
            norm_sumdecay=sumdecay/max(sumdecay)
            
            flux1=fftconvolve(pulse1, norm_sumdecay)
            flux2=fftconvolve(pulse2, norm_sumdecay)
        
            
            pulse3=fit_a3*np.exp(-.5*((time -fit_mu3)/fit_sig3)**2)
            decay3=np.exp(-(time)/fit_tau3)
            
            flux3=fftconvolve(pulse3, decay3)
            
            

            result[:,f-f_i]=results
            result_error[:,f-f_i]=results_error
            pulse1_spec[f-f_i,:]=pulse1[0:s_l]
            pulse2_spec[f-f_i,:]=pulse2[0:s_l]
            pulse3_spec[f-f_i,:]=pulse3[0:s_l]

            flux1_spec[f-f_i,:]=flux1[0:s_l]
            flux2_spec[f-f_i,:]=flux2[0:s_l]
            flux3_spec[f-f_i,:]=flux3[0:s_l]

            current_mdl= model(time,fit_a1,fit_mu1,fit_sig1,fit_tau1,fit_a2,fit_mu2,fit_sig2,fit_tau2,fit_a3,fit_mu3,fit_sig3,fit_tau3,fit_w)[0:s_l]
            diff=measured_signal2- current_mdl
            xi1= (diff)**2/(l*measured_signal2e)
            #xi1[56:60]=0
            xi=np.sum(xi1)    
            xirs[f-f_i]=xi/(s_l-len(ans))
               
            if plotter:
                fig,ax = plt.subplots(1,1,figsize=(16,10))

                ax.plot(pulse1+np.min(np.abs(measured_signal2)),'--',label="Fit Pulse 1",color='darkorange')
                ax.plot(pulse2+np.min(np.abs(measured_signal2)),'--', label="Fit Pulse 2", color='navy')
                ax.plot(pulse3+np.min(np.abs(measured_signal2)),'--', label="Fit Pulse 3", color='darkgreen')


                ax.plot(model(time,fit_a1,fit_mu1,fit_sig1,fit_tau1,0,fit_mu2,fit_sig2,fit_tau2,0,fit_mu3,fit_sig3,fit_tau3,fit_w), label="Fit Flux P1",color='darkorange')
                ax.plot(model(time,0,fit_mu1,fit_sig1,fit_tau1,fit_a2,fit_mu2,fit_sig2,fit_tau2,0,fit_mu3,fit_sig3,fit_tau3,fit_w), label="Fit Flux P2", color='navy')
                ax.plot(model(time,0,fit_mu1,fit_sig1,fit_tau1,0,fit_mu2,fit_sig2,fit_tau2,fit_a3,fit_mu3,fit_sig3,fit_tau3,fit_w), label="Fit Flux P3", color='darkgreen')

                ax.plot(measured_signal2,'.', label="Eovsa", color='red')
                ax.plot(current_mdl,color='black', label="Fit")
                if log_scale: plt.yscale("log")

                ax.set_xlabel('Time ('+t_unit+')'),plt.ylabel('Region Flux ('+ b_unit +')')
                ax.set_xlim(0,s_l)

                ax.set_title(str(fl[f])+ ''+ f_unit +'. Reduced Xi^2 = %.3f' %xirs[f-f_i])
                ax.legend(loc="upper right")
                plt.savefig(pref+'_3 pulse Sample_'+str(fl[f])+ ''+ f_unit +'.png')
                plt.close()
        except:
            result[:,f-f_i]= np.nan

    if plotter:

        try:

            t_w=(f_w/5) - 1
            ticks=np.round([t_w,2*t_w,3*t_w,4*t_w,5*t_w])
            ticks=ticks.astype('int32')

            fig,ax = plt.subplots(1,1,figsize=(16,10))
            ax.pcolormesh(pulse1_spec[0:f_w,:])
            ax.set_yticks(ticks)
            ax.set_yticklabels(fl[ticks+f_i])
            ax.set_title('Pulse 1')
            ax.set_ylabel('Freq ('+ f_unit +')')
            ax.set_xlabel('Time ('+t_unit+')')
            plt.savefig(pref+'_Pulse_1_spectrum_.png')
            plt.close()

            fig,ax = plt.subplots(1,1,figsize=(16,10))
            ax.pcolormesh(pulse2_spec[0:f_w,:])
            ax.set_yticks(ticks)
            ax.set_yticklabels(fl[ticks+f_i])
            ax.set_title('Pulse 2')
            ax.set_ylabel('Freq ('+ f_unit +')')
            ax.set_xlabel('Time ('+t_unit+')')
            plt.savefig(pref+'_Pulse_2_spectrum.png')
            plt.close()


            fig,ax = plt.subplots(1,1,figsize=(16,10))
            ax.pcolormesh(pulse3_spec[0:f_w,:])
            ax.set_yticks(ticks)
            ax.set_yticklabels(fl[ticks+f_i])
            ax.set_title('Pulse 3')
            ax.set_ylabel('Freq ('+ f_unit +')')
            ax.set_xlabel('Time ('+t_unit+')')
            plt.savefig(pref+'_Pulse_3_spectrum_.png')
            plt.close()


            fig,ax = plt.subplots(1,1,figsize=(16,10))
            ax.pcolormesh(flux1_spec[0:f_w,:])
            ax.set_yticks(ticks)
            ax.set_yticklabels(fl[ticks+f_i])
            ax.set_title('Pulse 1 Flux')
            ax.set_ylabel('Freq ('+ f_unit +')')
            ax.set_xlabel('Time ('+t_unit+')')
            plt.savefig(pref+'_Pulse_1_spectrumflux_.png')
            plt.close()

            fig,ax = plt.subplots(1,1,figsize=(16,10))
            ax.pcolormesh(flux2_spec[0:f_w,:])
            ax.set_yticks(ticks)
            ax.set_yticklabels(fl[ticks+f_i])
            ax.set_title('Pulse 2 Flux')
            ax.set_ylabel('Freq ('+ f_unit +')')
            ax.set_xlabel('Time ('+t_unit+')')
            plt.savefig(pref+'_Pulse_2_spectrumflux_.png')
            plt.close()

            fig,ax = plt.subplots(1,1,figsize=(16,10))
            ax.pcolormesh(flux3_spec[0:f_w,:])
            ax.set_yticks(ticks)
            ax.set_yticklabels(fl[ticks+f_i])
            ax.set_title('Pulse 3 Flux')
            ax.set_ylabel('Freq ('+ f_unit +')')
            ax.set_xlabel('Time ('+t_unit+')')
            plt.savefig(pref+'_Pulse_3_spectrumflux_.png')
            plt.close()

            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            ax.pcolormesh(pulse1_spec[0:f_w, :] + pulse2_spec[0:f_w, :]+ pulse3_spec[0:f_w, :])
            ax.set_yticks(ticks)
            ax.set_yticklabels(fl[ticks+f_i])
            ax.set_title('Pulse Composite')
            ax.set_ylabel('Freq ('+ f_unit +')')
            ax.set_xlabel('Time ('+t_unit+')')
            plt.savefig(pref + '_Pulse_spectrum_.png')
            plt.close()

            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            ax.pcolormesh(flux1_spec[0:f_w, :] + flux2_spec[0:f_w, :]+ flux3_spec[0:f_w, :])
            ax.set_yticks(ticks)
            ax.set_yticklabels(fl[ticks+f_i])
            ax.set_title('Pulse Flux Composite')
            ax.set_ylabel('Freq ('+ f_unit +')')
            ax.set_xlabel('Time ('+t_unit+')')
            plt.savefig(pref + '_Pulse_spectrumflux_.png')
            plt.close()

            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            ax.plot(xirs[0:f_w], label="Reduced Xi^2")
            ax.set_title('Reduced Xi^2')
            ax.legend(loc="upper right")
            ax.set_xticks(ticks)
            ax.set_xticklabels(fl[ticks+f_i])
            ax.set_xlabel('Freq ('+ f_unit +')')
            ax.set_ylabel('Reduced Xi^2')
            plt.savefig(pref + '_Reduced_XiSq_.png')
            plt.close()
            plt.clf()

            plt.rcParams['font.size'] = 12

            #### outlier scan #### (Own Function)

            lower_stop_a = min(lb[0], lb[4],lb[8])
            upper_stop_a = 1.05 * max(ub[0], ub[4], ub[8])

            lower_stop_mu = min(lb[1], lb[5],lb[9])
            upper_stop_mu = 1.05 * max(ub[1], ub[5], ub[9])

            lower_stop_sig = min(lb[2], lb[6],lb[10])
            upper_stop_sig = 1.05 * max(ub[2], ub[6], ub[10])

            lower_stop_tau = min(lb[3], lb[7],lb[11])
            upper_stop_tau = 1.05 * max(ub[3], ub[7], ub[11])

            lower_quartile_a = max(np.nanpercentile(result[0, :], 25), np.nanpercentile(result[4, :], 25),
                                   np.nanpercentile(result[8, :], 25))
            upper_quartile_a = max(np.nanpercentile(result[0, :], 75), np.nanpercentile(result[4, :], 75),
                                   np.nanpercentile(result[8, :], 75))

            lower_quartile_mu = max(np.nanpercentile(result[1, :], 25), np.nanpercentile(result[5, :], 25),
                                    np.nanpercentile(result[9, :], 25))
            upper_quartile_mu = max(np.nanpercentile(result[1, :], 75), np.nanpercentile(result[5, :], 75),
                                    np.nanpercentile(result[9, :], 75))

            lower_quartile_sig = max(np.nanpercentile(result[2, :], 25), np.nanpercentile(result[6, :], 25),
                                     np.nanpercentile(result[10, :], 25))
            upper_quartile_sig = max(np.nanpercentile(result[2, :], 75), np.nanpercentile(result[6, :], 75),
                                     np.nanpercentile(result[10, :], 75))

            lower_quartile_tau = max(np.nanpercentile(result[3, :], 25), np.nanpercentile(result[7, :], 25),
                                     np.nanpercentile(result[11, :], 25))
            upper_quartile_tau = max(np.nanpercentile(result[3, :], 75), np.nanpercentile(result[7, :], 75),
                                     np.nanpercentile(result[11, :], 75))

            avg_a_e = max(np.nanpercentile(result_error[0, :], 50), np.nanpercentile(result_error[4, :], 50),
                          np.nanpercentile(result_error[8, :], 50))

            avg_mu_e = max(np.nanpercentile(result_error[1, :], 50), np.nanpercentile(result_error[5, :], 50),
                           np.nanpercentile(result_error[9, :], 50))

            avg_sig_e = max(np.nanpercentile(result_error[2, :], 50), np.nanpercentile(result_error[6, :], 50),
                            np.nanpercentile(result_error[10, :], 50))

            avg_tau_e = max(np.nanpercentile(result_error[3, :], 50), np.nanpercentile(result_error[7, :], 50),
                            np.nanpercentile(result_error[11, :], 50))

            mid_a_max = max([np.nanpercentile(result[0, :], 50), np.nanpercentile(result[4, :], 50),
                             np.nanpercentile(result[8, :], 50)])
            mid_a_min = min([np.nanpercentile(result[0, :], 50), np.nanpercentile(result[4, :], 50),
                             np.nanpercentile(result[8, :], 50)])

            mid_mu_max = max([np.nanpercentile(result[1, :], 50), np.nanpercentile(result[5, :], 50),
                              np.nanpercentile(result[9, :], 50)])
            mid_mu_min = min([np.nanpercentile(result[1, :], 50), np.nanpercentile(result[5, :], 50),
                              np.nanpercentile(result[9, :], 50)])

            mid_sig_max = max([np.nanpercentile(result[2, :], 50), np.nanpercentile(result[6, :], 50),
                               np.nanpercentile(result[10, :], 50)])
            mid_sig_min = min([np.nanpercentile(result[2, :], 50), np.nanpercentile(result[6, :], 50),
                               np.nanpercentile(result[10, :], 50)])

            mid_tau_max = max([np.nanpercentile(result[3, :], 50), np.nanpercentile(result[7, :], 50),
                               np.nanpercentile(result[11, :], 50)])
            mid_tau_min = min([np.nanpercentile(result[3, :], 50), np.nanpercentile(result[7, :], 50),
                               np.nanpercentile(result[11, :], 50)])

            IQR_a = (upper_quartile_a - lower_quartile_a) * 2
            IQR_mu = (upper_quartile_mu - lower_quartile_mu) * 2
            IQR_sig = (upper_quartile_sig - lower_quartile_sig) * 2
            IQR_tau = (upper_quartile_tau - lower_quartile_tau) * 2

            a_min_1 = mid_a_min - IQR_a
            a_max_1 = mid_a_max + IQR_a

            mu_min_1 = mid_mu_min - IQR_mu
            mu_max_1 = mid_mu_max + IQR_mu

            sig_min_1 = mid_sig_min - IQR_sig
            sig_max_1 = mid_sig_max + IQR_sig

            tau_min_1 = mid_tau_min - IQR_tau
            tau_max_1 = mid_tau_max + IQR_tau

            a_min_e = mid_a_min - avg_a_e
            a_max_e = mid_a_max + avg_a_e

            mu_min_e = mid_mu_min - avg_mu_e
            mu_max_e = mid_mu_max + avg_mu_e

            sig_min_e = mid_sig_min - avg_sig_e
            sig_max_e = mid_sig_max + avg_sig_e

            tau_min_e = mid_tau_min - avg_tau_e
            tau_max_e = mid_tau_max + avg_tau_e

            a_min = max([lower_stop_a, a_min_1, a_min_e])
            a_max = min([upper_stop_a, a_max_1, a_max_e])

            mu_min = max([lower_stop_mu, mu_min_1, mu_min_e])
            mu_max = min([upper_stop_mu, mu_max_1, mu_max_e])

            sig_min = max([lower_stop_sig, sig_min_1, sig_min_e])
            sig_max = min([upper_stop_sig, sig_max_1, sig_max_e])

            tau_min = max([lower_stop_tau, tau_min_1, tau_min_e])
            tau_max = min([upper_stop_tau, tau_max_1, tau_max_e])
            ####

            fig, ax = plt.subplots(2, 2, figsize=(16, 10))


            ax[0][0].errorbar(fl[f_i:f_f], result[0, 0:f_w], yerr=result_error[0, 0:f_w], fmt='.', label="Amp 1",color='darkorange')
            ax[0][0].errorbar(fl[f_i:f_f], result[4, 0:f_w], yerr=result_error[4, 0:f_w], fmt='.', label="Amp 2", color='navy')
            ax[0][0].errorbar(fl[f_i:f_f],result[8,0:f_w],yerr=result_error[8,0:f_w], fmt='.',label="Amp 3", color='darkgreen')
            ax[0][0].legend(loc="upper right")
            ax[0][0].set_xticks(fl[ticks+f_i])
            ax[0][0].set_xticklabels(fl[ticks+f_i])
            ax[0][0].set_xlabel('Freq ('+ f_unit +')')
            ax[0][0].set_ylabel('Amplitude ('+ b_unit +'/s)')
            ax[0][0].set_ylim(a_min, a_max)

            ax[0][1].errorbar(fl[f_i:f_f], result[1, 0:f_w], yerr=result_error[1, 0:f_w], fmt='.', label="Mu 1",color='darkorange')
            ax[0][1].errorbar(fl[f_i:f_f], result[5, 0:f_w], yerr=result_error[5, 0:f_w], fmt='.', label="Mu 2", color='navy')
            ax[0][1].errorbar(fl[f_i:f_f],result[9,0:f_w],yerr=result_error[9,0:f_w], fmt='.',label="Mu 3", color='darkgreen')
            ax[0][1].legend(loc="upper left")
            ax[0][1].set_xticks(fl[ticks+f_i])
            ax[0][1].set_xticklabels(fl[ticks+f_i])
            ax[0][1].set_xlabel('Freq ('+ f_unit +')')
            ax[0][1].set_ylabel('Pulse Center ('+t_unit+')')
            ax[0][1].set_ylim(mu_min, mu_max)

            ax[1][0].errorbar(fl[f_i:f_f], result[2, 0:f_w], yerr=result_error[2, 0:f_w], fmt='.', label="Sigma 1",color='darkorange')
            ax[1][0].errorbar(fl[f_i:f_f], result[6, 0:f_w], yerr=result_error[6, 0:f_w], fmt='.', label="Sigma 2", color='navy')
            ax[1][0].errorbar(fl[f_i:f_f],result[10,0:f_w],yerr=result_error[10,0:f_w], fmt='.',label="Sigma 3", color='darkgreen')
            ax[1][0].legend(loc="upper left")
            ax[1][0].set_xticks(fl[ticks+f_i])
            ax[1][0].set_xticklabels(fl[ticks+f_i])
            ax[1][0].set_xlabel('Freq ('+ f_unit +')')
            ax[1][0].set_ylabel('Pulse Width ('+t_unit+')')
            ax[1][0].set_ylim(sig_min, sig_max)

            ax[1][1].errorbar(fl[f_i:f_f],result[3,0:f_w],yerr=result_error[3,0:f_w], fmt='.',label="Coupled Tau 1",color='darkorange')
            ax[1][1].errorbar(fl[f_i:f_f],result[7,0:f_w],yerr=result_error[7,0:f_w], fmt='.', label="Coupled Tau 2", color='navy')
            ax[1][1].errorbar(fl[f_i:f_f],result[11,0:f_w],yerr=result_error[11,0:f_w], fmt='.',label="Tau 3", color='darkgreen')
            ax[1][1].legend(loc="upper right")
            ax[1][1].set_xticks(fl[ticks+f_i])
            ax[1][1].set_xticklabels(fl[ticks+f_i])
            ax[1][1].set_xlabel('Freq ('+ f_unit +')')
            ax[1][1].set_ylabel('Decay Constant ('+t_unit+')')
            ax[1][1].set_ylim(tau_min, tau_max)

            plt.savefig(pref+'_pulse3_2mixedParameter_Panels_.png')
            plt.close()

            parms=[["freqs"],["fit_a1"], ["fit_mu1"], ["fit_sig1"], ["fit_tau1"],["fit_a2"], ["fit_mu2"], ["fit_sig2"], ["fit_tau2"], ["fit_a3"], ["fit_mu3"], ["fit_sig3"], ["fit_tau3"],["fit_w"],["e_fit_a1"], ["e_fit_mu1"], ["e_fit_sig1"], ["e_fit_tau1"],["e_fit_a2"], ["e_fit_mu2"], ["e_fit_sig2"], ["e_fit_tau2"], ["e_fit_a3"], ["e_fit_mu3"], ["e_fit_sig3"], ["e_fit_tau3"],["e_fit_w"],['reduced Xi^2']]
            result_2= np.append([fl[f_i:f_f]],result,0)
            result_2 = np.append(result_2, result_error,0)
            result_2 = np.append(result_2, [xirs], 0)
            result_2 = np.append(parms,result_2,1)
            np.savetxt(pref+"params.csv", result_2, delimiter=",", fmt="%s")

        except:
            print("Error, No fit found")

    else:
        return ans, xirs[1], cov