# This script is intended for PDN current profile generation, manipulation 
# It can support following type of waveforms.
### NOTE: CLK_Freq_in_GHz needs to inversely scaled by the factor used to scale Time_unit_in_sec 

#INFO: waveform_type A constant_clk:    Time_Length_in_ns | current_amplitude | current_floor < | clk conti| clk skip | this_freq_in_GHz >
#INFO: waveform_type B linear_slope_clk: Time_Length_in_ns | current_amplitude_start | current_amplitude_end |current_floor < | clk conti| clk skip | this_freq_in_GHz >
#INFO: waveform_type C scaled profile in .pkl:  "file path to envelope source (no space allowed)" | Time_unit_in_sec | Waveform_in_metric_unit | current_floor | time_scaling_factor | mag_scaling_factor | col_clk | col_data | skip_first_n_data  < | clk conti| clk skip | this_freq_in_GHz >
#INFO: waveform_type D scaled profile:  "file path to envelope source (no space allowed)" | Time_unit_in_sec | Waveform_in_metric_unit | current_floor | time_scaling_factor | mag_scaling_factor | skip_first_n_data  < | clk conti| clk skip | this_freq_in_GHz >
#INFO: waveform_type E random btw low/up bound: Time_Length_in_ns | current_low_bound | current_up_bound | current_floor  < | clk conti| clk skip | this_freq_in_GHz >
#INFO: waveform_type F linear_slope_no_clk: Time_Length_in_ns | current_amplitude_start | current_amplitude_end < | clk conti| clk skip | this_freq_in_GHz >

#INFO: CLK_Freq unit: GHz, it does NOT have to be logic freq, one can use lower number to reduce sampling freq
#INFO: 0. < CLK_DutyCycle < 1.
#INFO: 0. < CLK_T_RISE_as_ratio_of_CLK_Freq < 1.
#INFO: 0. < CLK_T_FALL_as_ratio_of_CLK_Freq < 1.
#INFO: 0. < (CLK_T_RISE_as_ratio_of_CLK_Freq + CLK_T_FALL_as_ratio_of_CLK_Freq) < CLK_DutyCycle

### START example input.params
##INFO: need to change VDD_in_Volt to 1.0, if waveform D is current vector, rather than power vector
# VDD_in_Volt  0.75
# PROFILE_POWER_0_OR_CURR_1   1
# CLK_Freq_in_GHz    3.0
# CLK_DutyCycle 0.9
# CLK_T_RISE_as_ratio_of_CLK_Freq 0.25
# CLK_T_FALL_as_ratio_of_CLK_Freq 0.25
# CLK_EDGE_EFF N

## usage 1
##INFO: Time_Length_in_ns  #Waveform_type  #Waveform_params
# A   10  5.  0.
# B   10  5.  100.    0.
# B   10  100. 55.    0.
# A   15 55. 0.
# C   C:\Users\jiangongwei\Documents\Python_data\NNE_power_trace_test1.pkl  1.e-9   1.     0.  1.0     0.25    Cycles   pNNE 0  1
# D   C:\Users\jiangongwei\Documents\Python_data\pwr_envelop.txt   1.e-9   1.0  0.  1.  1.  0  1
# A   10   55.0  0.
### NOTE: waveform C support .pkl to pick column titled "Cycles" and "pNNE", use 0.25 mag scaling

## usage 2 (clk gating)
##INFO: examples to activate clk gating (ccontinue 3 clks, skip 2 clks, while clk freq is 2GHz and 3.5GHz, respectively)
# B   10  0.  100.    0.    4   1   2.0
# B   10  100. 55.    0.    3   2   3.5
### END example input.params

import os, sys, getopt
import stat
from stat import S_IREAD, S_IRGRP, S_IROTH
# import math, cmath
# import shutil
import matplotlib.pyplot as plt
import scipy.interpolate
import numpy.random
import pandas as pd

file_dir = ""
file_in_para = ""
file_out_waveform = ""
voltage = 0.
waveform_params_list = []
is_print_wf = False
 
###
class CurrWaveform:
    currWaveform_list_time_ns = []  # unit: ns
    currWaveform_list_curr_Amp = []   # unit: Amp
    waveform_params_list = []

    profile_power_0_or_curr_1 = 0

    ### list of for scaling profiles, can support multiple profiles
    src_profile_envelope_fileName = ""
    src_profile_envelope_time_unit_in_sec = 1.
    src_profile_envelope_waveform_amp_unit = 1.

    clk_freq_norm = 0

    clk_freq = 0
    T_clk = 0
    T_clk_in_ns = 0

    clk_duty_cycle = 0.5
    clk_edge_peak_floor_ratio = 0.8
    t_rise_ratio_T = 0.1
    t_fall_ratio_T = 0.1
    nominal_I = 1.0
    I_lkg = 1.0
    voltage = 0.0
    r_curr_mag_scale_fac_charge_consv_h2 = 1.0
    r_curr_mag_scale_fac_charge_consv_h1 = r_curr_mag_scale_fac_charge_consv_h2 * 0.5  ### default
    
    is_clk_eff = True
    ### only applied to waveform C
    waveform_c_col_clk = ''
    waveform_c_col_data = ''
    waveform_c_time_scale_fac = 1.0 
    waveform_c_mag_scale_fac = 1.0 
    waveform_c_skip_n_data = 0
    #waveform_c_w_clk = True 
    ### only applied to waveform D
    waveform_d_time_scale_fac = 1.0 
    waveform_d_mag_scale_fac = 1.0 
    waveform_d_skip_n_data = 0
    #waveform_d_w_clk = True 

    ### these 3 paras are per waveform, will be overwritten by next waveform
    is_clk_gating = False
    numOfConsecutiveClk = 0
    numOfSkippedClk = 0

    currWaveform_list_time_ns.append(0)
    currWaveform_list_curr_Amp.append(0)

    ### Function: Read in waveform parameters
    def __init__(self, file_in_para):
        
        line_cnt = 0
        wf_cnt = 0
        set_cnt = 0
        with open(r'%s'%file_in_para, 'r') as fin:
            cln_str = fin.readline()
            while cln_str:
                cln_str_src = cln_str.lstrip(' ').rstrip('\n')
                if cln_str_src == '' or cln_str_src[0] == '#':
                    cln_str = fin.readline()
                    continue 

                cln_str = cln_str_src.split()
                if cln_str[0] == 'CLK_Freq_in_GHz':
                    self.clk_freq = float( cln_str[1]) * 1.e9 
                    self.clk_freq_norm = self.clk_freq
                    self.T_clk = 1. / self.clk_freq
                    self.T_clk_in_ns = self.T_clk * 1.e9
                    print('#INFO: CLK freq (GHz):\t' + str(self.clk_freq / 1.e9))
                    print('#INFO: CLK Cycle (ns):\t' + str(self.T_clk_in_ns))
                    set_cnt = set_cnt +1
                    
                elif cln_str[0] == 'CLK_DutyCycle':
                    self.clk_duty_cycle = float( cln_str[1])
                    print('#INFO: CLK duty cycle:\t' + str(self.clk_duty_cycle))
                    set_cnt = set_cnt +1

                elif cln_str[0] == 'CLK_T_RISE_as_ratio_of_CLK_Freq':
                    self.t_rise_ratio_T = float( cln_str[1])
                    print('#INFO: CLK rise time (ratio of T):\t' + str(self.t_rise_ratio_T))
                    set_cnt = set_cnt +1
                    
                elif cln_str[0] == 'CLK_T_FALL_as_ratio_of_CLK_Freq':
                    self.t_fall_ratio_T = float( cln_str[1])
                    print('#INFO: CLK fall time (ratio of T):\t' + str(self.t_rise_ratio_T))
                    set_cnt = set_cnt +1
                    
                elif cln_str[0] == 'VDD_in_Volt':
                    self.voltage = float(cln_str[1])
                    if abs(self.voltage) < 1.e-6:
                        print('\n#ERRORR: Vdd rail voltage needs to be larger than 0, input value is ' + str(self.voltage) + '\n\n')
                        exit(-1)
                    print('#INFO: Vdd rail voltage (V):\t' + str(self.voltage))
                    set_cnt = set_cnt +1
                
                elif cln_str[0] == 'PROFILE_POWER_0_OR_CURR_1':
                    self.profile_power_0_or_curr_1 = int(cln_str[1])
                    set_cnt = set_cnt +1

                elif cln_str[0] == 'CLK_EDGE_EFF':
                    if cln_str[1] == 'Y' or cln_str[1] == 'y':
                        self.is_clk_eff = True 
                    elif cln_str[1] == 'N' or cln_str[1] == 'n':
                        self.is_clk_eff = False
                    else:
                        print('#ERROR: CLK EDGE EFF definition has to be Y/y or N/n, quit. \n')
                        exit(-1)
                    set_cnt = set_cnt +1
                    print('#INFO: CLK edge effect:', cln_str[1])

                else:  
                    ## only proceed if general settings are completed
                    if set_cnt != 7:
                        print('#ERROR: 7 general settings needed, %d identified, quit\n', set_cnt)
                        exit(-1)
                    # read in waveform params
                    self.waveform_params_list.append(cln_str_src)
                    wf_cnt = wf_cnt + 1
                    print('#INFO: Reading waveform def. ' + str(wf_cnt) + ': ' + self.waveform_params_list[-1])
                
                cln_str = fin.readline()
        fin.close()

        ### sanity check
        if (self.t_rise_ratio_T + self.t_fall_ratio_T) > self.clk_duty_cycle:
            print ('#ERROR: summation of clk rise time and fall time ratio cannot exceed duty cycle ratio\n')
            sys.exit(1)
        self.r_curr_mag_scale_fac_charge_consv_h2 =  1./ (self.clk_duty_cycle - 0.5 * self.t_rise_ratio_T - 0.5 * self.t_fall_ratio_T)  ### this is assuming floor is 0
        print('\n#INFO: defaullt waveforms WITH clk edges with be scaled by ' + str(self.r_curr_mag_scale_fac_charge_consv_h2) + '\n')

    def ClkGatingClear(self):
        self.is_clk_gating = False
        self.numOfConsecutiveClk = 0
        self.numOfSkippedClk = 0

    def ClearWaveformInfo(self):
        self.src_profile_envelope_fileName = ''
        self.src_profile_envelope_time_unit_in_sec = 1.
        self.src_profile_envelope_waveform_amp_unit = 1.
        I_floor_ratio = 0.
        self.waveform_d_time_scale_fac = 1.
        self.waveform_d_mag_scale_fac = 1.

        self.waveform_c_col_clk = ''
        self.waveform_c_col_data = ''
        self.waveform_c_time_scale_fac = 1.0 
        self.waveform_c_mag_scale_fac = 1.0 
        self.waveform_c_skip_n_data = 0
        #self.waveform_c_w_clk = True 
        ### only applied to waveform D
        self.waveform_d_time_scale_fac = 1.0 
        self.waveform_d_mag_scale_fac = 1.0 
        self.waveform_d_skip_n_data = 0
        #self.waveform_d_w_clk = True 

    def ReadClkGatingInfo(self, numOfConsecutiveClk_, numOfSkippedClk_, waveFormName):
        self.is_clk_gating = True 
        self.numOfConsecutiveClk = numOfConsecutiveClk_
        self.numOfSkippedClk = numOfSkippedClk_
        print("#INFO: waveform type "+ waveFormName+" clk gating enabled, consecutive/skipped clk = " + str(self.numOfConsecutiveClk) + ' / ' + str(self.numOfSkippedClk))

    def UpdateFreq(self, freq_in_GHz):
        self.clk_freq = freq_in_GHz * 1.e9
        self.T_clk = 1./self.clk_freq
        self.T_clk_in_ns = self.T_clk * 1.e9

    def RestoreFreq(self):
        self.clk_freq = self.clk_freq_norm
        self.T_clk = 1./self.clk_freq
        self.T_clk_in_ns = self.T_clk * 1.e9      

    def AddDelayWF(self, delay_wf_length_ns, delay_wf_curr):
        tStart = self.currWaveform_list_time_ns[-1]
        self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns * self.t_rise_ratio_T)
        self.currWaveform_list_curr_Amp.append(delay_wf_curr)

        self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns * self.t_rise_ratio_T + delay_wf_length_ns)
        self.currWaveform_list_curr_Amp.append(delay_wf_curr)

    ### Function: Add one unit
    ### Optional clk gating: parameters after "I_floor_ratio" is optional to enable clk gating
    def AddOneUnit(self, I_amp, I_floor_val, clk_seq = 0):
        tStart = self.currWaveform_list_time_ns[-1]
        ### no clk gating
        if not self.is_clk_gating:
            self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns * self.t_rise_ratio_T)
            self.currWaveform_list_curr_Amp.append(I_amp)

            self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns * self.clk_duty_cycle - self.T_clk_in_ns*self.t_fall_ratio_T)
            self.currWaveform_list_curr_Amp.append(I_amp)

            self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns * self.clk_duty_cycle)
            self.currWaveform_list_curr_Amp.append(I_floor_val)

            self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns)
            self.currWaveform_list_curr_Amp.append(I_floor_val)
            return ### exit function if no clk gating

        ### continue if clk gated 
        rem = clk_seq % (self.numOfConsecutiveClk + self.numOfSkippedClk)
        if rem < self.numOfConsecutiveClk:
            tStart = self.currWaveform_list_time_ns[-1]
            self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns * self.t_rise_ratio_T)
            self.currWaveform_list_curr_Amp.append(I_amp)

            self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns * self.clk_duty_cycle - self.T_clk_in_ns*self.t_fall_ratio_T)
            self.currWaveform_list_curr_Amp.append(I_amp)

            self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns * self.clk_duty_cycle)
            self.currWaveform_list_curr_Amp.append(I_floor_val)

            self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns)
            self.currWaveform_list_curr_Amp.append(I_floor_val)
        else:
            self.AddOneUnit(I_floor_val, I_floor_val)   ### ??? gated clk pwr to be determined 

    def UpdateClkEdgeWaveformPara(self, I_floor_ratio):
        self.r_curr_mag_scale_fac_charge_consv_h2 =  1./ ( I_floor_ratio + (self.clk_duty_cycle - 0.5 * self.t_rise_ratio_T - 0.5 * self.t_fall_ratio_T) * (1. - I_floor_ratio) )  ### this is assuming floor is 0
        self.r_curr_mag_scale_fac_charge_consv_h1 = self.r_curr_mag_scale_fac_charge_consv_h2 * I_floor_ratio
        print('#INFO: this waveforms with be scaled (if clk edged) by: ' + str(self.r_curr_mag_scale_fac_charge_consv_h2) + ' and ' + str(self.r_curr_mag_scale_fac_charge_consv_h1) +'\n')

    ### Function: Append time length of nominal clk current
    def AddConstCLK(self, NumOfUnit, I_amp):
        for i in range (0, NumOfUnit):
            if self.is_clk_eff:
                self.AddOneUnit(I_amp * self.r_curr_mag_scale_fac_charge_consv_h2, I_amp * self.r_curr_mag_scale_fac_charge_consv_h1, i)
            else: ## if no clk eff, this is really just piece wise linear
                self.AddOneUnit(I_amp, I_amp, i)

    ### Function: Append linear ramp up curent from I_start to I_end witihn t_ramp time
    def AddLinearSlopeCurr(self, numOfUnit, I_start, I_end):
        I_step = (I_end - I_start)/numOfUnit
        for i in range(0, numOfUnit):
            I_tmp = I_start + (i+1) * I_step
            if self.is_clk_eff:
                self.AddOneUnit(I_tmp * self.r_curr_mag_scale_fac_charge_consv_h2, I_tmp * self.r_curr_mag_scale_fac_charge_consv_h1, i)
            else:
                self.AddOneUnit(I_tmp, I_tmp, i)
    
    def AddLinearSlopeCurr_noClk(self, numOfUnit, I_start, I_end):
        ### Note: no clk is involved, hence no current scaling 
        I_step = (I_end - I_start)/numOfUnit
        for i in range(0, numOfUnit):
            I_tmp = I_start + (i+1) * I_step
            tStart = self.currWaveform_list_time_ns[-1]
            self.currWaveform_list_time_ns.append(tStart + 0.25 * self.T_clk_in_ns)
            self.currWaveform_list_curr_Amp.append(I_tmp + 0.25 * I_step)
            self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns)
            self.currWaveform_list_curr_Amp.append(I_tmp + I_step)           

    def AddScalingCurr_waveformC(self):
        ### read in source file
        if not os.path.exists(self.src_profile_envelope_fileName):
            print('#ERROR: Source profile envelope file <' + self.src_profile_envelope_fileName + '> does not exist !')
            sys.exit(1)

        if self.profile_power_0_or_curr_1 == 0:
            print('#INFO: Input is power value, waveform C magnitdue is scaled by 1/voltage = '+ str(1/self.voltage))
        if self.waveform_c_skip_n_data != 0:
            print('#INFO: skipping first ' + str(self.waveform_c_skip_n_data) + ' data samplings')

        df = pd.read_pickle(self.src_profile_envelope_fileName)
        src_profile_time_in_ns = df[self.waveform_c_col_clk] 
        src_profile_time_in_ns = src_profile_time_in_ns[self.waveform_c_skip_n_data : -1] - src_profile_time_in_ns[self.waveform_c_skip_n_data]  ## time starts a 0

        ### NOTE:  obsoleted: use Cycles to scale and obtain time steps
        #src_profile_time_in_ns = src_profile_time_in_ns * self.src_profile_envelope_time_unit_in_sec * 1.e9 * self.T_clk_in_ns * self.waveform_c_time_scale_fac

        # get time step directly from pkl
        src_profile_time_in_ns = src_profile_time_in_ns * self.src_profile_envelope_time_unit_in_sec * 1.e9

        src_profile_amplitude = df[self.waveform_c_col_data] 
        src_profile_amplitude = src_profile_amplitude[self.waveform_c_skip_n_data : -1]
        if self.profile_power_0_or_curr_1 == 0:
            src_profile_amplitude = src_profile_amplitude * self.src_profile_envelope_waveform_amp_unit/ self.voltage * self.waveform_c_mag_scale_fac
        else:
            src_profile_amplitude = src_profile_amplitude * self.src_profile_envelope_waveform_amp_unit * self.waveform_c_mag_scale_fac
        print('#INFO: waveform C columns: \n')
        print(df.columns)
        #tmp_col = df.columns
        #tmp_row0= df.iloc[0,:]
        src_profile_time_in_ns = src_profile_time_in_ns.values.tolist()
        src_profile_amplitude  = src_profile_amplitude.values.tolist()
        func_intep = scipy.interpolate.interp1d(src_profile_time_in_ns, src_profile_amplitude)
        curr_delta_min = 1.e8
        curr_delta_max = -1.e-8
        amp_last = 0; 
        numOfUnit = int( (src_profile_time_in_ns[-1] - src_profile_time_in_ns[0]) / self.T_clk_in_ns )
        for i in range (0, numOfUnit):
            time_ = i * self.T_clk_in_ns
            amp_ = func_intep(time_)

            ##### temp skip some units
            # skip_adj = self.src_profile_envelope_time_unit_in_sec * 1.e12
            # # print('DeBuG: skip_adj = ', skip_adj)
            # if time_ > 5000. * skip_adj and  time_ < 24000 * skip_adj:
            #     continue
            # elif time_ > 34000. * skip_adj and  time_ < 40000 * skip_adj:
            #     continue
            # elif time_> 60000 * skip_adj:
            #     continue 
            ##### temp skip some units end 

            #### current amplitude scaled such that charge is const
            if self.is_clk_eff :
                self.AddOneUnit(amp_ * self.r_curr_mag_scale_fac_charge_consv_h2, amp_ * self.r_curr_mag_scale_fac_charge_consv_h1, i)  
            else:
                self.AddOneUnit(amp_, amp_, i)  
            #self.AddOneUnit( src_profile_amplitude[i] * self.r_curr_mag_scale_fac_charge_consv_h2, I_floor_ratio, i) 
            curr_delta = amp_ - amp_last
            if curr_delta > curr_delta_max:
                curr_delta_max = curr_delta
            if curr_delta < curr_delta_min:
                curr_delta_min = curr_delta
            amp_last = amp_ 

        print('#INFO: scaled worst pos. delta current per cycle(before charge scaling): '+str(curr_delta_max) + ', scaled worst neg. delta current per cycle(before charge scaling): ' + str(curr_delta_min) )

    def AddScalingCurr_waveformD(self):
        ### read in source file
        if not os.path.exists(self.src_profile_envelope_fileName):
            print('#ERROR: Source profile envelope file <' + self.src_profile_envelope_fileName + '> does not exist !')
            sys.exit(1)

        if self.profile_power_0_or_curr_1 == 0:
            print('#INFO: Input is power value, waveform D magnitdue is scaled by 1/voltage = '+ str(1/self.voltage))
        if self.waveform_d_skip_n_data != 0:
            print('#INFO: skipping first ' + str(self.waveform_d_skip_n_data) + ' data samplings')

        src_profile_time_in_ns = []
        src_profile_amplitude  = []

        with open(r'%s'%self.src_profile_envelope_fileName, 'r') as fin:
            cln_str = fin.readline()
            time_ST = 0
            time_st_fnd = False
            data_cnt = 0
            while cln_str:
                cln_str_src = cln_str.lstrip(' ').rstrip('\n')
                if cln_str_src == '' or cln_str_src[0] == '#' or 'time' in cln_str_src or 'TIME' in cln_str_src or 'Time' in cln_str_src: ### skip empty line or lines start with #, or 'time'
                    cln_str = fin.readline()
                    continue    
                if ',' in cln_str:
                    cln_str = cln_str_src.split(',') 
                else:   
                    cln_str = cln_str_src.split()
                
                time_ns = float(cln_str[0]) * self.src_profile_envelope_time_unit_in_sec * 1.e9
                data_cnt = data_cnt + 1

                if data_cnt >= self.waveform_d_skip_n_data: 
                    if not time_st_fnd:
                        time_ST = time_ns
                        time_st_fnd = True

                    src_profile_time_in_ns.append( self.waveform_d_time_scale_fac * ( time_ns - time_ST))    ### nominal profile starts from 0
                    ### Note: divided by voltage to obtian current
                    if self.profile_power_0_or_curr_1 == 0:
                        curr_ = float(cln_str[1]) * self.src_profile_envelope_waveform_amp_unit/ self.voltage
                    else:
                        curr_ = float(cln_str[1]) * self.src_profile_envelope_waveform_amp_unit

                    src_profile_amplitude.append( curr_ * self.waveform_d_mag_scale_fac )     

                cln_str = fin.readline()                

        fin.close()
        #### interpolate 
        time_len = src_profile_time_in_ns[-1] - src_profile_time_in_ns[0]
        numOfUnit = int( time_len / self.T_clk_in_ns)
        func_intep = scipy.interpolate.interp1d(src_profile_time_in_ns, src_profile_amplitude)

        curr_delta_min = 1.e8
        curr_delta_max = -1.e-8
        amp_last = 0; 

        for i in range (0, numOfUnit):
            time_ = i * self.T_clk_in_ns
            amp_ = func_intep(time_)

            ##### temp skip some units
            # skip_adj = self.src_profile_envelope_time_unit_in_sec * 1.e12
            skip_adj = 1
            # print('DeBuG: skip_adj = ', skip_adj)
            
            if time_ > 1000. * skip_adj and  time_ < 24000 * skip_adj:
                continue
            elif time_ > 29000. * skip_adj and  time_ < 32000 * skip_adj:
                continue
            elif time_ > 35000. * skip_adj and  time_ < 38000 * skip_adj:
                continue
            elif time_> 48000 * skip_adj:
                continue 
            ##### temp skip some units end 


            #### current amplitude scaled such that charge is const
            if self.is_clk_eff:
                self.AddOneUnit(amp_ * self.r_curr_mag_scale_fac_charge_consv_h2, amp_ * self.r_curr_mag_scale_fac_charge_consv_h1, i)      
            else:
                self.AddOneUnit(amp_, amp_, i)  

            curr_delta = amp_ - amp_last
            if curr_delta > curr_delta_max:
                curr_delta_max = curr_delta
            if curr_delta < curr_delta_min:
                curr_delta_min = curr_delta
            amp_last = amp_ 

        print('#INFO: scaled worst pos. delta current per cycle(before charge scaling): '+str(curr_delta_max) + ',scaled worst neg. delta current per cycle(before charge scaling): ' + str(curr_delta_min) )

    ### Function: Add random current within given upper and lower bound
    def AddRandWithinRange(self, numOfUnit, I_bd_lo, I_bd_up):
        rng = numpy.random.default_rng()
        currValList = rng.random((numOfUnit))
        for idx, i in enumerate( currValList):
            curr = I_bd_lo + (I_bd_up - I_bd_lo) * i
            if self.is_clk_eff:
                self.AddOneUnit(curr * self.r_curr_mag_scale_fac_charge_consv_h2, curr * self.r_curr_mag_scale_fac_charge_consv_h1, idx)
            else:
                self.AddOneUnit(curr, curr, idx)

    ### Function: compose the actual waveform based on parameters
    def CompositeWaveform(self):
        for wfp in self.waveform_params_list:
            wfp_orig = wfp
            wfp = wfp.split() ### split by space
            time_wf_ns = 0.

            self.ClkGatingClear()

            if wfp[0] != 'D' and wfp[0] != 'C' :
                time_wf_ns = float( wfp[1])
            numOfUnit =  int(time_wf_ns/ self.T_clk_in_ns)
        
            if wfp[0] == 'A':
                if len(wfp) < 4:
                    print("#ERROR: waveform type A parameters insufficient: " + wfp_orig)
                    sys.exit(1)
                else:
                    I_curr = float( wfp[2])
                    I_floor_ratio = float( wfp[3])

                    if len(wfp) == 7: ## clk gated 
                        self.ReadClkGatingInfo(int( wfp[4]), int( wfp[5]), 'A')
                        self.UpdateFreq(float(wfp[6]))
                    self.UpdateClkEdgeWaveformPara(I_floor_ratio)
                    self.AddConstCLK(numOfUnit, I_curr )
                self.ClkGatingClear()
                self.RestoreFreq()

            elif wfp[0] == 'B':
                if len(wfp) < 5:
                    print("#ERROR: waveform type B parameters insufficient: " + wfp_orig)
                    sys.exit(1)
                else:
                    I_start = float( wfp[2])
                    I_end = float( wfp[3])
                    I_floor_ratio = float( wfp[4])          

                    if len(wfp) == 8:
                        self.ReadClkGatingInfo(int( wfp[5]), int( wfp[6]), 'B')
                        self.UpdateFreq(float(wfp[7]))

                    self.UpdateClkEdgeWaveformPara(I_floor_ratio)
                    self.AddLinearSlopeCurr(numOfUnit, I_start, I_end)               
                self.ClkGatingClear()
                self.RestoreFreq()

            elif wfp[0] == 'C': 
                if len(wfp) < 10:
                    print("#ERROR: waveform type C parameters insufficient: " + wfp_orig)
                    sys.exit(1)
                else:
                    self.src_profile_envelope_fileName = wfp[1]
                    self.src_profile_envelope_time_unit_in_sec = float(wfp[2]) 
                    self.src_profile_envelope_waveform_amp_unit = float(wfp[3])
                    I_floor_ratio = float(wfp[4])
                    self.waveform_c_time_scale_fac = float(wfp[5])
                    self.waveform_c_mag_scale_fac  = float(wfp[6])
                    ### column names in pkl
                    self.waveform_c_col_clk = wfp[7]
                    self.waveform_c_col_data = wfp[8]
                    self.waveform_c_skip_n_data = int(wfp[9])

                    if len(wfp) == 13:
                        self.ReadClkGatingInfo(int( wfp[10]), int( wfp[11]), 'D')
                        self.UpdateFreq(float(wfp[12]))

                    self.UpdateClkEdgeWaveformPara(I_floor_ratio)
                    self.AddScalingCurr_waveformC()
                self.ClkGatingClear()
                self.ClearWaveformInfo()
                self.RestoreFreq()
                
            elif wfp[0] == 'D': 
                if len(wfp) < 8:
                    print("#ERROR: waveform type D parameters insufficient: " + wfp_orig)
                    sys.exit(1)
                else:
                    self.src_profile_envelope_fileName = wfp[1]
                    self.src_profile_envelope_time_unit_in_sec = float(wfp[2]) 
                    self.src_profile_envelope_waveform_amp_unit = float(wfp[3])
                    I_floor_ratio = float(wfp[4])
                    self.waveform_d_time_scale_fac = float(wfp[5])
                    self.waveform_d_mag_scale_fac  = float(wfp[6])
                    self.waveform_d_skip_n_data = int(wfp[7])

                    if len(wfp) == 11:
                        self.ReadClkGatingInfo(int( wfp[8]), int( wfp[9]), 'D')
                        self.UpdateFreq(float(wfp[10]))

                    self.UpdateClkEdgeWaveformPara(I_floor_ratio)
                    self.AddScalingCurr_waveformD()
                self.ClkGatingClear()
                self.ClearWaveformInfo()
                self.RestoreFreq()

            elif wfp[0] == 'E':
                if len(wfp) < 5:
                    print('#ERROR: waveform type E parameters insufficient: ' + wfp_orig)
                    sys.exit(1)
                else:
                    I_floor_ratio = float(wfp[4])
                    I_bd_lo = float(wfp[2])
                    I_bd_up = float(wfp[3])
                    if I_bd_lo > I_bd_up:
                        tmp = I_bd_lo
                        I_bd_lo = I_bd_up
                        I_bd_up = tmp 

                    if len(wfp) == 8:
                        self.ReadClkGatingInfo(int( wfp[5]), int( wfp[6]), 'E')
                        self.UpdateFreq(float(wfp[7]))

                    self.UpdateClkEdgeWaveformPara(I_floor_ratio)
                    self.AddRandWithinRange(numOfUnit, I_bd_lo, I_bd_up)
                self.ClkGatingClear()
                self.RestoreFreq()

            elif wfp[0] == 'F':
                if len(wfp) < 4:
                    print("#ERROR: waveform type F parameters insufficient: " + wfp_orig)
                    sys.exit(1)
                else:
                    I_start = float( wfp[2])
                    I_end = float( wfp[3])        

                    if len(wfp) == 7:
                        self.ReadClkGatingInfo(int( wfp[4]), int( wfp[5]), 'F')
                        self.UpdateFreq(float(wfp[6]))

                    self.AddLinearSlopeCurr_noClk(numOfUnit, I_start, I_end)
                self.ClkGatingClear()
                self.RestoreFreq()
            elif wfp[0] == 'G':
                if len(wfp) < 2:
                    print("#ERROR: waveform type G parameters insufficient: " + wfp_orig)
                    sys.exit(1)
                else:
                    delay_wf_length_ns = float( wfp[1])
                    delay_wf_curr = float( wfp[2])
                    self.AddDelayWF(delay_wf_length_ns, delay_wf_curr )


    def WriteWaveform(self, fileName):
        fout = open(fileName, 'w+')
        #fout.write('*Time(s) \t Current(Amp)\n')
        #fout.write('I_src_???  YourNode1??? YourNode2??? pwl \n')

        fout.write('.subckt    curr_src_???? \n')
        fout.write('+ pin_pos ref_gnd\n\n')
        fout.write('I_currSrc    pin_pos    ref_gnd    pwl\n')

        leng_rcd = len( self.currWaveform_list_time_ns)
        for i in range(0, leng_rcd):
            fout.write('+ ' + str(self.currWaveform_list_time_ns[i]) + 'e-9\t' + str(self.currWaveform_list_curr_Amp[i]) + '\t\n')  
        fout.write('\n.ends\n')
        fout.close()
        print('#INFO: Current waveform output as PWL to ' + fileName)

    def WriteWaveform_InTimFormat(self, fileName):
        fileName_ = fileName.split('.')
        if fileName_[-1] != 'tim' or fileName_[-1] != 'TIM': 
            fileName = fileName + '.tim'

        fout = open(fileName, 'w+')
        leng_rcd = len( self.currWaveform_list_time_ns)
        fout.write('BEGIN  TIMEDATA' + '\n')
        fout.write('# T ( NSEC  V  R 50 )' + '\n')
        fout.write('%   time      voltage' + '\n')  ### Note: this is not typo, "voltage" is just format for TIM file
        for i in range(0, leng_rcd):
            fout.write(str(self.currWaveform_list_time_ns[i]) + '\t' + str(self.currWaveform_list_curr_Amp[i]) + '\n')  
        fout.write('END' + '\n')
        fout.close()
        print('#INFO: Current waveform TIM format output: ' + fileName)

    def WriteWaveform_ToCSV(self, fileName):
        fileName = fileName + '.csv'
        fout = open(fileName, 'w+')
        leng_rcd = len( self.currWaveform_list_time_ns)
        for i in range(0, leng_rcd):
            fout.write(str(self.currWaveform_list_time_ns[i]) + 'e-9, \t' + str(self.currWaveform_list_curr_Amp[i]) + '\t\n')  

        fout.close()
        print('#INFO: Current waveform output as PWL to ' + fileName)

    def WriteWaveform_ToSimplis(self, fileName):
        fileName = fileName + '_simplis.csv'

        if os.path.isfile(fileName):
            os.chmod(fileName, stat.S_IWRITE)
            os.remove(fileName)

        fout = open(fileName, 'w+')
        leng_rcd = len( self.currWaveform_list_time_ns)
        if leng_rcd > 1e5:
            print('\n\n#WARNING: SIMPLIS cannot take more than 1e5 data points ! ' + str(leng_rcd) +' data points exported \n\n')

        fout.write('START_DATA SHIFT_FIRST_TO_ZERO FORMAT=CSV,\n')
        for i in range(0, leng_rcd):
            fout.write(str(self.currWaveform_list_time_ns[i]) + 'e-9, "\t' + str(self.currWaveform_list_curr_Amp[i]) + ' "\t\n')  

        fout.close()
        # NOTE: simplis output fomrat files are made read only due to SIMPLIS tends to change file in run
        
        os.chmod(fileName, S_IREAD|S_IRGRP|S_IROTH)
        print('#INFO: Current waveform output as PWL to READ ONLY file : ' + fileName)

    def PlotWaveform(self):
        plt.rcParams.update({'font.size': 15})
        plt.plot(self.currWaveform_list_time_ns, self.currWaveform_list_curr_Amp, 'b', linewidth=0.2, label='Waveform')
        plt.legend(loc=1)
        plt.xlabel('Time (ns)', fontsize=14, weight='bold')
        plt.ylabel('Current(Amp)', fontsize=14, weight='bold')
        plt.title('Current Profle')

        ax = plt.subplot(111)
        plt.grid(True, which='both')	
        ax.grid(which='major', linewidth=0.8)
        ax.grid(which='minor', linestyle=':', linewidth=0.5)
        ax.minorticks_on()        
        plt.savefig(file_dir + 'current_waveform.png', dpi = 600)
        plt.show()
        plt.clf()
        plt.close()
	

# Main function 
try:
	opts,args = getopt.getopt(sys.argv[1:],'d:i:o:p')
except getopt.GetoptError:
	print('\nUsage: python pdn_current_gene_main.py [-d file directory] [-i input.params] [-o out_curr_profile.tim] [-p <if waveforms print>]')
	sys.exit(2)
if (not opts) and args:
	print('\nUsage: python pdn_current_gene_main.py [-d file directory] [-i input.params] [-o out_curr_profile.tim] [-p <if waveforms print>]')
	sys.exit(2)

for o,a in opts:
    if o =='-d':
        if os.name == 'nt':
            file_dir = a.lstrip(' ').rstrip(' ') + '\\'
        else:
            file_dir = a.lstrip(' ').rstrip(' ') + '/'
    if o =='-i':
        file_in_para = a.lstrip(' ').rstrip(' ')
    if o =='-o':
        file_out_waveform = a.lstrip(' ').rstrip(' ')
    if o == '-p':
        is_print_wf = True

# BEGIN read input parameters
file_in_para = file_dir + file_in_para
file_out_waveform = file_dir + file_out_waveform

if os.path.exists(file_in_para):
    print('#INFO: input parameters file:    ' + file_in_para)
    print('#INFO: output waveform TIM file: ' + file_out_waveform)
else:
    print('#ERROR: intput parameter file <' + file_in_para + '> does not exist !')
    sys.exit(1)


# END read input parameters
print('\n#INFO: Start to run waveform generation !')

waveformInst = CurrWaveform(file_in_para)
waveformInst.CompositeWaveform()
# waveformInst.WriteWaveform(file_out_waveform)
# waveformInst.WriteWaveform_InTimFormat(file_out_waveform)
waveformInst.WriteWaveform_ToCSV(file_out_waveform)
waveformInst.WriteWaveform_ToSimplis(file_out_waveform)

if is_print_wf:
    waveformInst.PlotWaveform()

print ("#INFO: Waveform generation exit normally")