

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import pulses1 as p1
import pulses1_m as p1_m
import pulses2_separate as p2_s
import pulses2_2mixed as p2_2m
import pulses3_separate as p3_s
import pulses3_2mixed as p3_2m
import pulses3_3mixed as p3_3m
import pulses4_separate as p4_s
import pulses4_2mixed as p4_2m
import pulses4_2_2mixed as p4_2_2m
import pulses4_3mixed as p4_3m
import pulses4_4mixed as p4_4m

def model_1pulse(t,a1,mu1,sig1,tau1):

    
    pulse1=a1*np.exp(-.5*((t-mu1)/sig1)**2)
    decay1=np.exp(-t / tau1)
    flux1=np.convolve(pulse1, decay1)
    

    test_signal= flux1

    t_l = len(t)

    test_signal = flux1[0:t_l]
    
    return test_signal





###############################################################################################################

def sorter(m_s, m_s_e, freqs, p_times,guess, min_g, max_g, t_i, t_f, f):

    s_l= t_f - t_i


    xirs=np.zeros([13])
    xirs[0:13]=np.nan

    model_results = np.zeros([12,18])
    model_results[:,:] = np.nan

    model_flag=np.zeros([12])
    model_flag[0:12]=0.0






    init_guess_1pulse = guess[0:4]
    init_guess_1pulse_mixed = np.append(guess[0:4], [guess[7], guess[16]])
    init_guess_2pulse_separate = guess[0:8]
    init_guess_2pulse_mixed = np.append(guess[0:8],1)
    init_guess_3pulse_separate= guess[0:12]
    init_guess_3pulse_2mixed= np.append(guess[0:12],1)
    init_guess_3pulse_3mixed = np.append(guess[0:11],1)
    init_guess_4pulse_separate = guess[0:16]
    init_guess_4pulse_2mixed= guess[0:17]
    init_guess_4pulse_2_2mixed= guess[0:18]
    init_guess_4pulse_3mixed = np.append(guess[0:15],guess[16])
    init_guess_4pulse_4mixed = np.append(guess[0:11],[guess[12],guess[13],guess[14],guess[16]])

    lb_1pulse = min_g[0:4]
    lb_1pulse_mixed = np.append(min_g[0:4], [min_g[7], min_g[16]])
    lb_2pulse_separate = min_g[0:8]
    lb_2pulse_mixed = np.append(min_g[0:8],min_g[16])
    lb_3pulse_separate=min_g[0:12]
    lb_3pulse_2mixed= np.append(min_g[0:12],min_g[16])
    lb_3pulse_3mixed = np.append(min_g[0:11],min_g[16])
    lb_4pulse_separate = min_g[0:16]
    lb_4pulse_2mixed= np.append(min_g[0:16],min_g[16])
    lb_4pulse_2_2mixed= np.append(min_g[0:16],[min_g[16],min_g[17]])
    lb_4pulse_3mixed = np.append(min_g[0:15],min_g[16])
    lb_4pulse_4mixed = np.append(min_g[0:11],[min_g[12],min_g[13],min_g[14],min_g[16]])

    ub_1pulse = max_g[0:4]
    ub_1pulse_mixed = np.append(max_g[0:4], [max_g[7], max_g[16]])
    ub_2pulse_separate = max_g[0:8]
    ub_2pulse_mixed = np.append(max_g[0:8],max_g[16])
    ub_3pulse_separate=max_g[0:12]
    ub_3pulse_2mixed= np.append(max_g[0:12],max_g[16])
    ub_3pulse_3mixed = np.append(max_g[0:11],max_g[16])
    ub_4pulse_separate = max_g[0:16]
    ub_4pulse_2mixed= np.append(max_g[0:16],max_g[16])
    ub_4pulse_2_2mixed= np.append(max_g[0:16],[max_g[16],max_g[17]])
    ub_4pulse_3mixed = np.append(max_g[0:15],max_g[16])
    ub_4pulse_4mixed = np.append(max_g[0:11],[max_g[12],max_g[13],max_g[14],max_g[16]])


    oddflag= False
    l=1






    measured_signal2 = m_s[t_i:t_f, f]
    measured_signal2e = m_s_e[t_i:t_f, f]

    temp_m_s =np.zeros([s_l,3])
    temp_m_s_e =np.zeros([s_l,3])
    temp_m_s[:,0]  = m_s[t_i:t_f, f]
    temp_m_s_e[:,0] = m_s_e[t_i:t_f, f]
    temp_m_s[:,1]  = m_s[t_i:t_f, f]
    temp_m_s_e[:,1] = m_s_e[t_i:t_f, f]
    temp_m_s[:,2]  = m_s[t_i:t_f, f]
    temp_m_s_e[:,2] = m_s_e[t_i:t_f, f]


    time= p_times[t_i:t_f]-p_times[t_i]

    if s_l % 2 != 0:
        s_l = s_l + 1
        oddflag=True



    if s_l != np.size(temp_m_s,0):

        t1 =np.array([[m_s[t_f, f]],[m_s[t_f, f]],[m_s[t_f, f]]])
        te =np.array([[m_s_e[t_f, f]],[m_s_e[t_f, f]],[m_s_e[t_f, f]]])
        temp_m_s = np.append(temp_m_s, np.transpose(t1),axis=0)
        temp_m_s_e = np.append(temp_m_s_e, np.transpose(te),axis=0)
        time= np.append(time,time[s_l-2])



    try:

        ans_1pulse, xirs[0] = p1.p1(m_s=temp_m_s, m_s_e=temp_m_s_e, freqs=freqs, p_times=p_times, f_i=0, f_f=2, t_i=0, t_f=s_l,
              ig=init_guess_1pulse, ub=ub_1pulse, lb=lb_1pulse, plotter=False)
        model_results[0, 0:(len(ans_1pulse))] = ans_1pulse



    except:
        xirs[0]=np.nan
        model_flag[0]=np.nan



    try:
        ans_1pulse_mixed, xirs[11]= p1_m.p1_m(m_s=temp_m_s, m_s_e=temp_m_s_e, freqs=freqs, p_times=p_times, f_i=0, f_f=2, t_i=0, t_f=s_l,
                  ig=init_guess_1pulse_mixed, ub=ub_1pulse_mixed, lb=lb_1pulse_mixed,plotter=False)
        model_results[11,0:(len(ans_1pulse_mixed))]=ans_1pulse_mixed

    except:
        xirs[11]=np.nan
        model_flag[11]=np.nan



    try:

        ans_2pulse_separate,xirs[1] = p2_s.p2_s(m_s=temp_m_s, m_s_e=temp_m_s_e, freqs=freqs,p_times=p_times, f_i=0, f_f=2, t_i=0, t_f=s_l, ig=init_guess_2pulse_separate, ub=ub_2pulse_separate, lb=lb_2pulse_separate, plotter=False)
        model_results[1, 0:(len(ans_2pulse_separate))] = ans_2pulse_separate

    except:
        xirs[1]=np.nan
        model_flag[1]=np.nan


    try:

        ans_2pulse_mixed,xirs[2] = p2_2m.p2_2m(m_s=temp_m_s,m_s_e=temp_m_s_e,freqs=freqs,p_times=p_times,f_i=0,f_f=2,t_i=0, t_f=s_l,ig=init_guess_2pulse_mixed,ub=ub_2pulse_mixed,lb=lb_2pulse_mixed,plotter=False)
        model_results[2, 0:(len(ans_2pulse_mixed))] = ans_2pulse_mixed

    except:
        xirs[2]=np.nan
        model_flag[2]=np.nan



    try:
        ans_3pulse_separate,xirs[3] = p3_s.p3_s(m_s=temp_m_s,m_s_e=temp_m_s_e,freqs=freqs,p_times=p_times,f_i=0,f_f=2,t_i=0, t_f=s_l,ig=init_guess_3pulse_separate,ub=ub_3pulse_separate,lb=lb_3pulse_separate,plotter=False)
        model_results[3, 0:(len(ans_3pulse_separate))] = ans_3pulse_separate

    except:
        xirs[3]=np.nan
        model_flag[3]=np.nan


    try:
        ans_3pulse_2mixed,xirs[4]= p3_2m.p3_2m(m_s=temp_m_s,m_s_e=temp_m_s_e,freqs=freqs,p_times=p_times,f_i=0,f_f=2,t_i=0, t_f=s_l,ig=init_guess_3pulse_2mixed,ub=ub_3pulse_2mixed,lb=lb_3pulse_2mixed,plotter=False)
        model_results[4, :(len(ans_3pulse_2mixed))] = ans_3pulse_2mixed


    except:
        xirs[4]=np.nan
        model_flag[4]=np.nan


    try:
        ans_3pulse_3mixed,xirs[5] = p3_3m.p3_3m(m_s=temp_m_s,m_s_e=temp_m_s_e,freqs=freqs,p_times=p_times,f_i=0,f_f=2,t_i=0, t_f=s_l,ig=init_guess_3pulse_3mixed,ub=ub_3pulse_3mixed,lb=lb_3pulse_3mixed,plotter=False)
        model_results[5, 0:(len(ans_3pulse_3mixed))] = ans_3pulse_3mixed

    except:
        xirs[5]=np.nan
        model_flag[5]=np.nan


    try:

        ans_4pulse_separate, xirs[6] = p4_s.p4_s(m_s=temp_m_s,m_s_e=temp_m_s_e,freqs=freqs,p_times=p_times,f_i=0,f_f=2,t_i=0, t_f=s_l,ig=init_guess_4pulse_separate,ub=ub_4pulse_separate,lb=lb_4pulse_separate,plotter=False)
        model_results[6, 0:(len(ans_4pulse_separate))] = ans_4pulse_separate

    except:
        xirs[6]=np.nan
        model_flag[6]=np.nan


    try:
        ans_4pulse_2mixed, xirs[7] = p4_2m.p4_2m(m_s=temp_m_s,m_s_e=temp_m_s_e,freqs=freqs,p_times=p_times,f_i=0,f_f=2,t_i=0, t_f=s_l,ig=init_guess_4pulse_2mixed,ub=ub_4pulse_2mixed,lb=lb_4pulse_2mixed,plotter=False)
        model_results[7, :(len(ans_4pulse_2mixed))] = ans_4pulse_2mixed

    except:
        xirs[7]=np.nan
        model_flag[7]=np.nan

    l=1
    try:

        ans_4pulse_2_2mixed,xirs[8] = p4_2_2m.p4_2_2m(m_s=temp_m_s,m_s_e=temp_m_s_e,freqs=freqs,p_times=p_times,f_i=0,f_f=2,t_i=0, t_f=s_l,ig=init_guess_4pulse_2_2mixed,ub=ub_4pulse_2_2mixed,lb=lb_4pulse_2_2mixed,plotter=False)
        model_results[8, 0:(len(ans_4pulse_2_2mixed))] = ans_4pulse_2_2mixed


    except:
        xirs[8]=np.nan
        model_flag[8]=np.nan



    try:

        ans_4pulse_3mixed,xirs[9] = p4_3m.p4_3m(m_s=temp_m_s,m_s_e=temp_m_s_e,freqs=freqs,p_times=p_times,f_i=0,f_f=2,t_i=0, t_f=s_l,ig=init_guess_4pulse_3mixed,ub=ub_4pulse_3mixed,lb=lb_4pulse_3mixed,plotter=False)
        model_results[9, 0:(len(ans_4pulse_3mixed))] = ans_4pulse_3mixed

    except:
        xirs[9]=np.nan
        model_flag[9]=np.nan



    try:

        ans_4pulse_4mixed,xirs[10] = p4_4m.p4_4m(m_s=temp_m_s,m_s_e=temp_m_s_e,freqs=freqs,p_times=p_times,f_i=0,f_f=2,t_i=0, t_f=s_l,ig=init_guess_4pulse_4mixed,ub=ub_4pulse_4mixed,lb=lb_4pulse_4mixed,plotter=False)
        model_results[10, :(len(ans_4pulse_4mixed))] = ans_4pulse_4mixed


    except:
        xirs[10]=np.nan
        model_flag[10]=np.nan



    good_models = np.where(np.logical_and(xirs >= np.nanmin(xirs[np.nonzero(xirs)]), xirs <= 1.2 * np.nanmin(xirs[np.nonzero(xirs)])))
    best_model = np.where(xirs == np.nanmin(xirs[np.nonzero(xirs)]))
    fit_model=[]
    if best_model==[]: best_model=[0]
    try:
        best_model=best_model[0]
    except:
        pass

    try:
        best_model[0][0]
        flag=1
    except:
        flag=0



    l=1
    if flag==0:
        if best_model[0]==[0]:
            fit_model= ans_1pulse
        if best_model[0]==[11]:
            fit_model= ans_1pulse_mixed
        if best_model[0]==[1]:
            fit_model= ans_2pulse_separate
        if best_model[0]==[2]:
            fit_model= ans_2pulse_mixed
        if best_model[0]==[3]:
            fit_model= ans_3pulse_separate
        if best_model[0]==[4]:
            fit_model= ans_3pulse_2mixed
        if best_model[0]==[5]:
            fit_model= ans_3pulse_3mixed
        if best_model[0]==[6]:
            fit_model= ans_4pulse_separate
        if best_model[0]==[7]:
            fit_model= ans_4pulse_2mixed
        if best_model[0]==[8]:
            fit_model= ans_4pulse_2_2mixed
        if best_model[0]==[9]:
            fit_model= ans_4pulse_3mixed
        if best_model[0]==[10]:
            fit_model= ans_4pulse_4mixed
    else:
        pass

    if flag==1:
        if best_model[0][:]==[0]:
            fit_model= ans_1pulse
        if best_model[0][:]==[11]:
            fit_model= ans_1pulse_mixed
        if best_model[0][:]==[1]:
            fit_model= ans_2pulse_separate
        if best_model[0][:]==[2]:
            fit_model= ans_2pulse_mixed
        if best_model[0][:]==[3]:
            fit_model= ans_3pulse_separate
        if best_model[0][:]==[4]:
            fit_model= ans_3pulse_2mixed
        if best_model[0][:]==[5]:
            fit_model= ans_3pulse_3mixed
        if best_model[0][:]==[6]:
            fit_model= ans_4pulse_separate
        if best_model[0][:]==[7]:
            fit_model= ans_4pulse_2mixed
        if best_model[0][:]==[8]:
            fit_model= ans_4pulse_2_2mixed
        if best_model[0][:]==[9]:
            fit_model= ans_4pulse_3mixed
        if best_model[0][:]==[10]:
            fit_model= ans_4pulse_4mixed
    else:
        pass







    if np.any(np.isin(good_models,0)):
        print("1 Pulse is the good fit with Reduced Xi^2 of ",xirs[0])
        model_flag[0]=1

    if np.any(np.isin(good_models,11)):
        print("1 Pulse with 2 decays is the good fit with Reduced Xi^2 of ",xirs[11])
        model_flag[11]=1

    if np.any(np.isin(good_models,1)):
        print("2 Separate Pulses is the good fit with Reduced Xi^2 of ",xirs[1])
        model_flag[1]=1

    if np.any(np.isin(good_models,2)):
        print("2 Pulses with mixed decays  is the good fit with Reduced Xi^2 of ",xirs[2])
        model_flag[2]=1

    if np.any(np.isin(good_models,3)):
        print("3 Separate Pulses is the good fit with Reduced Xi^2 of ",xirs[3])
        model_flag[3]=1

    if np.any(np.isin(good_models,4)):
        print("3 Pulses 2 of which have mixed decays is the good fit with Reduced Xi^2 of ",xirs[4])
        model_flag[4]=1

    if np.any(np.isin(good_models,5)):
        print("3 Pulses all of which have mixed decays is the good fit with Reduced Xi^2 of ",xirs[5])
        model_flag[5]=1

    if np.any(np.isin(good_models,6)):
        print("4 Separate Pulses is the good fit with Reduced Xi^2 of ",xirs[6])
        model_flag[6]=1

    if np.any(np.isin(good_models,7)):
        print("4 Pulses 2 of which have mixed decays is the good fit with Reduced Xi^2 of ",xirs[7])
        model_flag[7]=1

    if np.any(np.isin(good_models,8)):
        print("4 Pulses with 2 sets of 2 mixed decays is the good fit with Reduced Xi^2 of ",xirs[8])
        model_flag[8]=1

    if np.any(np.isin(good_models,9)):
        print("4 Pulses 3 of which have mixed decays is the good fit with Reduced Xi^2 of ",xirs[9])
        model_flag[9]=1

    if np.any(np.isin(good_models,10)):
        print("4 Pulses all of which have mixed decays is the good fit with Reduced Xi^2 of ",xirs[10])
        model_flag[10]=1


    plt.rcParams['font.size'] = 12

    return best_model, xirs, fit_model, model_results


