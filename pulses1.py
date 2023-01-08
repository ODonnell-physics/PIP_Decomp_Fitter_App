# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 11:47:56 2022

@author: Brian
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def model(t, a1, mu1, sig1, tau1):
    try:
        measured_signal2
    except:
        measured_signal2 = [0, 0]
    pulse1 = a1 * np.exp(-.5 * ((t - mu1) / sig1) ** 2)
    decay1 = np.exp(-t / tau1)
    t_l = len(t)
    flux1 = np.convolve(pulse1, decay1,mode='full')

    test_signal = flux1[0:t_l] + np.min(np.abs(measured_signal2))

    return test_signal


def p1(m_s, m_s_e, freqs, p_times, f_i, f_f, t_i, t_f, ig, ub, lb, plotter=True, f_unit="GHz", b_unit="SFU", t_unit="s",log_scale=False):


    f_w = f_f - f_i
    s_l = t_f - t_i
    s_l_=s_l
    if s_l % 2 != 0 : s_l_=s_l+1

    pulse1_spec = np.zeros([f_w, s_l_])
    flux1_spec = np.zeros([f_w, s_l_])

    def model(t, a1, mu1, sig1, tau1):

        pulse1 = a1 * np.exp(-.5 * ((t - mu1) / sig1) ** 2)
        decay1 = np.exp(-t / tau1)
        flux1 = np.convolve(pulse1, decay1,mode='full')
        t_l=len(t)

        test_signal = flux1[0:t_l] + np.min(np.abs(measured_signal2))

        return test_signal

    fl = np.around(freqs, 2)

    if plotter:
        if not os.path.exists('./1_pulse_fit'):
            os.makedirs('./1_pulse_fit')

    pref = './1_pulse_fit/1_'

    xirs = np.zeros([f_w])
    result = np.zeros([4, f_w])
    result_error = np.zeros([4, f_w])
    oddflag= False




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

            ans, cov = fit
            fit_a1, fit_mu1, fit_sig1, fit_tau1 = ans
            fit_sa1, fit_smu1, fit_ssig1, fit_stau1 = np.sqrt(np.diag(cov))

            results = np.array([fit_a1, fit_mu1, fit_sig1, fit_tau1])
            results_error = np.array([fit_sa1, fit_smu1, fit_ssig1, fit_stau1])

            result[:, f - f_i] = results
            result_error[:, f - f_i] = results_error

            pulse1 = fit_a1 * np.exp(-.5 * ((time - fit_mu1) / fit_sig1) ** 2)
            flux1 = model(time, fit_a1, fit_mu1, fit_sig1, fit_tau1)

            pulse1_spec[f-f_i, 0:s_l] = pulse1[0:s_l]
            flux1_spec[f-f_i, 0:s_l] = flux1[0:s_l]

            current_mdl = model(time, fit_a1, fit_mu1, fit_sig1, fit_tau1)[0:s_l]
            diff = measured_signal2[0:s_l] - current_mdl
            xi1 = diff ** 2 / (l * measured_signal2e[0:s_l])
            # xi1[56:60]=0
            xir = np.sum(xi1) / (s_l - len(ans))
            xirs[f-f_i] = xir



            if plotter:
                plt.plot(pulse1 + np.min(np.abs(measured_signal2)),'--', label="Fit Pulse 1", color='darkorange')
                plt.plot(measured_signal2,'.', label="Eovsa", color='red')
                plt.plot(flux1, color='black', label="Fit")
                plt.xlabel('Time ('+t_unit+')'), plt.ylabel('Region Flux ('+ b_unit +')')
                plt.xlim(0, s_l)
                if log_scale: plt.yscale("log")
                plt.title(str(fl[f]) + ''+ f_unit +' reduced Xi^2 = %.3f' % xir)
                plt.legend(loc="upper right")
                plt.savefig(pref + str(freqs[f]) + ''+ f_unit +'.png')
                plt.close()
                plt.cla()
                plt.clf()
        except:
            pass
        #### outlier scan ####

    if plotter:
        if oddflag:
            s_l=s_l-1
        lower_stop_a = lb[0]
        upper_stop_a = 1.05 * ub[0]

        lower_stop_mu = lb[1]
        upper_stop_mu = 1.05 * ub[1]

        lower_stop_sig = lb[2]
        upper_stop_sig = 1.05 * ub[2]

        lower_stop_tau = lb[3]
        upper_stop_tau = 1.05 * ub[3]

        lower_quartile_a = np.nanpercentile(result[0, :], 25)
        upper_quartile_a = np.nanpercentile(result[0, :], 75)

        lower_quartile_mu = np.nanpercentile(result[1, :], 25)
        upper_quartile_mu = np.nanpercentile(result[1, :], 75)

        lower_quartile_sig = np.nanpercentile(result[2, :], 25)
        upper_quartile_sig = np.nanpercentile(result[2, :], 75)

        lower_quartile_tau = np.nanpercentile(result[3, :], 25)
        upper_quartile_tau = np.nanpercentile(result[3, :], 75)



        mid_a_max = np.nanpercentile(result[0, :], 50)
        mid_a_min = np.nanpercentile(result[0, :], 50)

        mid_mu_max = np.nanpercentile(result[1, :], 50)
        mid_mu_min = np.nanpercentile(result[1, :], 50)
        mid_sig_max = np.nanpercentile(result[2, :], 50)
        mid_sig_min = np.nanpercentile(result[2, :], 50)

        mid_tau_max = np.nanpercentile(result[3, :], 50)
        mid_tau_min = np.nanpercentile(result[3, :], 50)

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



        a_min = max([lower_stop_a, a_min_1])
        a_max = min([upper_stop_a, a_max_1])

        mu_min = max([lower_stop_mu, mu_min_1])
        mu_max = min([upper_stop_mu, mu_max_1])

        sig_min = max([lower_stop_sig, sig_min_1])
        sig_max = min([upper_stop_sig, sig_max_1])

        tau_min = max([lower_stop_tau, tau_min_1])
        tau_max = min([upper_stop_tau, tau_max_1])




        try:
            t_w = (f_w / 5) - 1
            ticks = np.round([t_w, 2 * t_w, 3 * t_w, 4 * t_w, 5 * t_w])
            ticks = ticks.astype('int32')

            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            ax.pcolormesh(pulse1_spec[0:f_w, 0:s_l])
            ax.set_yticks(ticks)
            ax.set_yticklabels(fl[ticks+f_i])
            ax.set_title('Pulse 1')
            ax.set_ylabel('Freq ('+ f_unit +')')
            ax.set_xlabel('Time ('+t_unit+')')
            plt.savefig(pref + '_Pulse_1_spectrum.png')
            plt.close()

            fig, ax = plt.subplots(1, 1, figsize=(16, 10))
            ax.pcolormesh(flux1_spec[0:f_w, 0:s_l])
            ax.set_yticks(ticks)
            ax.set_yticklabels(fl[ticks+f_i])
            ax.set_title('Pulse 1 Flux')
            ax.set_ylabel('Freq ('+ f_unit +')')
            ax.set_xlabel('Time ('+t_unit+')')
            plt.savefig(pref + '_Pulse_1_spectrumflux_.png')
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

            xr = np.arange(0, f_w)

            fig, ax = plt.subplots(2, 2, figsize=(16, 10))

            ax[0][0].errorbar(fl[f_i:f_f],result[0, 0:f_w], yerr=result_error[0, 0:f_w], fmt='.', label="Amp", color='darkorange')
            ax[0][0].legend(loc="upper right")
            ax[0][0].set_xticks(fl[ticks+f_i])
            ax[0][0].set_xticklabels(fl[ticks+f_i])
            ax[0][0].set_xlabel('Freq ('+ f_unit +')')
            ax[0][0].set_ylabel('Amplitude ('+ b_unit +'/s)')
            ax[0][0].set_ylim(a_min, a_max)

            ax[0][1].errorbar(fl[f_i:f_f], result[1, 0:f_w], yerr=result_error[1, 0:f_w], fmt='.', label="Mu", color='darkorange')
            ax[0][1].legend(loc="upper left")
            ax[0][1].set_xticks(fl[ticks+f_i])
            ax[0][1].set_xticklabels(fl[ticks+f_i])
            ax[0][1].set_xlabel('Freq ('+ f_unit +')')
            ax[0][1].set_ylabel('Pulse Center ('+t_unit+')')
            ax[0][1].set_ylim(mu_min, mu_max)

            ax[1][0].errorbar(fl[f_i:f_f], result[2, 0:f_w], yerr=result_error[2, 0:f_w], fmt='.', label="Sigma", color='darkorange')
            ax[1][0].legend(loc="upper left")
            ax[1][0].set_xticks(fl[ticks+f_i])
            ax[1][0].set_xticklabels(fl[ticks+f_i])
            ax[1][0].set_xlabel('Freq ('+ f_unit +')')
            ax[1][0].set_ylabel('Pulse Width ('+t_unit+')')
            ax[1][0].set_ylim(sig_min, sig_max)

            ax[1][1].errorbar(fl[f_i:f_f], result[3, 0:f_w], yerr=result_error[3, 0:f_w], fmt='.', label="Tau", color='darkorange')
            ax[1][1].legend(loc="upper right")
            ax[1][1].set_xticks(fl[ticks+f_i])
            ax[1][1].set_xticklabels(fl[ticks+f_i])
            ax[1][1].set_xlabel('Freq ('+ f_unit +')')
            ax[1][1].set_ylabel('Decay Constant ('+t_unit+')')
            ax[1][1].set_ylim(tau_min, tau_max)

            plt.savefig(pref + '_1_pulse_Parameter_Panels_.png')
            plt.close()
            plt.clf()

            parms=[["freqs"],["fit_a1"], ["fit_mu1"], ["fit_sig1"], ["fit_tau1"],["e_fit_a1"], ["e_fit_mu1"], ["e_fit_sig1"], ["e_fit_tau1"],['reduced Xi^2']]
            result_2= np.append([fl[f_i:f_f]],result,0)
            result_2 = np.append(result_2, result_error,0)
            result_2 = np.append(result_2, [xirs], 0)
            result_2 = np.append(parms,result_2,1)
            np.savetxt(pref+"params.csv", result_2, delimiter=",", fmt="%s")


        except:
            print(" ")

    else:
        return ans, xir

