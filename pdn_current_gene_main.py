# This script is intended for PDN current profile generation, manipulation 
# It can support following type of waveforms.
#INFO: waveform_type A constant_clk:    Time_Length_in_ns | current_amplitude | current_floor
#INFO: waveform_type B linear_slope:    Time_Length_in_ns | current_amplitude_start | current_amplitude_end |current_floor
#INFO: waveform_type C clock_gating:    Time_Length_in_ns | current_amplitude | current_floor | Num_of_consecutive_clks | Num_of_skipped_clks
#INFO: waveform_type D scaled profile:  "file path to envelope source (no space allowed)" | Time_unit_in_sec | Waveform_unit_if_convert_to_unit_1
#INFO: 0. < CLK_DutyCycle < 1.
#INFO: 0. < CLK_T_RISE_as_ratio_of_CLK_Freq < 1.
#INFO: 0. < CLK_T_FALL_as_ratio_of_CLK_Freq < 1.
#INFO: 0. < (CLK_T_RISE_as_ratio_of_CLK_Freq + CLK_T_FALL_as_ratio_of_CLK_Freq) < CLK_DutyCycle

### START example input.params
# VDD_in_Volt  0.75
# CLK_Freq_in_GHz    3.0
# CLK_DutyCycle 0.5
# CLK_T_RISE_as_ratio_of_CLK_Freq 0.07
# CLK_T_FALL_as_ratio_of_CLK_Freq 0.07

# #INFO: Time_Length_in_ns  #Waveform_type  #Waveform_params
# A   10  5.  0.
# B   10  5.  100.    0.
# B   10  100. 55.    0.
# A   15 55. 0.
# C   20 60. 5.  2   3
# D   C:\Users\jiangongwei\Documents\Python_data\pwr_envelop.txt   1.e-9   1.0
# A   10   55.0  0.
### END example input.params

import os, sys, getopt
# import math, cmath
# import shutil
import matplotlib.pyplot as plt
import scipy.interpolate
import numpy.random

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

    ### list of for scaling profiles, can support multiple profiles
    src_profile_envelope_fileName = ""
    src_profile_envelope_time_unit_in_sec = 1.
    src_profile_envelope_waveform_unit = 1.

    clk_freq = 0
    T_clk = 0
    T_clk_in_ns = 0
    clk_duty_cycle = 0.5
    t_rise_ratio_T = 0.1
    t_fall_ratio_T = 0.1
    nominal_I = 1.0
    I_lkg = 1.0
    voltage = 0.0

    currWaveform_list_time_ns.append(0)
    currWaveform_list_curr_Amp.append(0)

    ### Function: Read in waveform parameters
    def __init__(self, file_in_para):
        
        line_cnt = 0
        wf_cnt = 0
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
                    self.T_clk = 1. / self.clk_freq
                    self.T_clk_in_ns = self.T_clk * 1.e9
                    print('#INFO: CLK freq (GHz):\t' + str(self.clk_freq / 1.e9))
                    print('#INFO: CLK Cycle (ns):\t' + str(self.T_clk_in_ns))
                    
                elif cln_str[0] == 'CLK_DutyCycle':
                    self.clk_duty_cycle = float( cln_str[1])
                    print('#INFO: CLK duty cycle:\t' + str(self.clk_duty_cycle))
                    
                elif cln_str[0] == 'CLK_T_RISE_as_ratio_of_CLK_Freq':
                    self.t_rise_ratio_T = float( cln_str[1])
                    print('#INFO: CLK rise time (ratio of T):\t' + str(self.t_rise_ratio_T))
                    
                elif cln_str[0] == 'CLK_T_FALL_as_ratio_of_CLK_Freq':
                    self.t_fall_ratio_T = float( cln_str[1])
                    print('#INFO: CLK fall time (ratio of T):\t' + str(self.t_rise_ratio_T))
                    
                elif cln_str[0] == 'VDD_in_Volt':
                    self.voltage = float(cln_str[1])
                    print('#INFO: Vdd rail voltage (V):\t' + str(self.voltage))

                else:  # read in waveform params
                    self.waveform_params_list.append(cln_str_src)
                    wf_cnt = wf_cnt + 1
                    print('#INFO: Reading waveform def. ' + str(wf_cnt) + ': ' + self.waveform_params_list[-1])
                
                cln_str = fin.readline()
                ### sanity check
                if (self.t_rise_ratio_T + self.t_fall_ratio_T) > self.clk_duty_cycle:
                    print ('#ERROR: summation of clk rise time and fall time ratio cannot exceed duty cycle ratio\n')
                    sys.exit(1)
        fin.close()


    ### Function: Add one unit
    def AddOneUnit(self, I_amp, I_floor):
        tStart = self.currWaveform_list_time_ns[-1]
        self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns * self.t_rise_ratio_T)
        self.currWaveform_list_curr_Amp.append(I_amp)

        self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns * self.clk_duty_cycle - self.T_clk_in_ns*self.t_fall_ratio_T)
        self.currWaveform_list_curr_Amp.append(I_amp)

        self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns * self.clk_duty_cycle)
        self.currWaveform_list_curr_Amp.append(I_floor)

        self.currWaveform_list_time_ns.append(tStart + self.T_clk_in_ns)
        self.currWaveform_list_curr_Amp.append(I_floor)


    ### Function: Append time length of nominal clk current
    def AddConstCLK(self, NumOfUnit, I_amp, I_floor):
        for i in range (0, NumOfUnit):
            self.AddOneUnit(I_amp, I_floor )

    ### Function: Append linear ramp up curent from I_start to I_end witihn t_ramp time
    def AddLinearSlopeCurr(self, numOfUnit, I_start, I_end, I_floor):
        I_step = (I_end - I_start)/numOfUnit
        for i in range(0, numOfUnit):
            I_tmp = I_start + (i+1) * I_step
            self.AddOneUnit(I_tmp, I_floor)


    ### Function: Append current demonstrated gated clk cycles, could be multiple patterns
    ### skip numOfSkippedClk, after numOfConsecutiveClk
    def AddClkGatingCurr(self, numOfUnit, I_amp, I_floor, numOfConsecutiveClk, numOfSkippedClk):
        for i in range (0, numOfUnit):
            rem = i % (numOfConsecutiveClk + numOfSkippedClk)
            if rem < numOfConsecutiveClk:
                self.AddOneUnit(I_amp, I_floor)
            else:
                self.AddOneUnit(I_floor, I_floor)

    def AddScalingCurr(self):
        ### read in source file
        if not os.path.exists(self.src_profile_envelope_fileName):
            print('#ERROR: Source profile envelope file <' + self.src_profile_envelope_fileName + '> does not exist !')
            sys.exit(1)

        src_profile_time_in_ns = []
        src_profile_amplitude  = []

        with open(r'%s'%self.src_profile_envelope_fileName, 'r') as fin:
            cln_str = fin.readline()
            time_ST = 0
            time_st_fnd = False
            while cln_str:
                cln_str_src = cln_str.lstrip(' ').rstrip('\n')
                if cln_str_src == '' or cln_str_src[0] == '#':
                    cln_str = fin.readline()
                    continue    

                cln_str = cln_str_src.split()   
                time_ns = float(cln_str[0]) * self.src_profile_envelope_time_unit_in_sec * 1.e9
                if not time_st_fnd:
                    time_ST = time_ns
                    time_st_fnd = True

                src_profile_time_in_ns.append(time_ns - time_ST)    ### nominal profile starts from 0
                ### Note: divided by voltage to obtian current
                curr_ = float(cln_str[1]) * self.src_profile_envelope_waveform_unit/ self.voltage
                src_profile_amplitude.append( curr_ )     
                cln_str = fin.readline()
        fin.close()
        #### interpolate 
        time_len = src_profile_time_in_ns[-1] - src_profile_time_in_ns[0]
        numOfUnit = int( time_len / self.T_clk_in_ns)
        func_intep = scipy.interpolate.interp1d(src_profile_time_in_ns, src_profile_amplitude)

        for i in range (0, numOfUnit):
            time_ = i * self.T_clk_in_ns
            amp_ = func_intep(time_)
            I_floor = 0.
            self.AddOneUnit(amp_, I_floor)      

    ### Function: Add random current within given upper and lower bound
    def AddRandWithinRange(self, numOfUnit, I_floor, I_bd_lo, I_bd_up):
        rng = numpy.random.default_rng()
        currValList = rng.random((numOfUnit,))
        for i in currValList:
            self.AddOneUnit(I_bd_lo + (I_bd_up - I_bd_lo) * i, I_floor)

    ### Function: compose the actual waveform based on parameters
    def CompositeWaveform(self):
        for wfp in self.waveform_params_list:
            wfp = wfp.split()
            time_wf_ns = 0.
            if wfp[0] != 'D':
                time_wf_ns = float( wfp[1])
            numOfUnit =  int(time_wf_ns/ self.T_clk_in_ns)
        
            if wfp[0] == 'A':
                if len(wfp) < 4:
                    print("#ERROR: waveform type A parameters insufficient: " + wfp)
                    sys.exit(1)
                else:
                    I_curr = float( wfp[2])
                    I_floor = float( wfp[3])
                    self.AddConstCLK(numOfUnit, I_curr, I_floor)
            elif wfp[0] == 'B':
                if len(wfp) < 5:
                    print("#ERROR: waveform type B parameters insufficient: " + wfp)
                    sys.exit(1)
                else:
                    I_start = float( wfp[2])
                    I_end = float( wfp[3])
                    I_floor = float( wfp[4])          
                    self.AddLinearSlopeCurr(numOfUnit, I_start, I_end, I_floor)
            elif wfp[0] == 'C':
                if len(wfp) < 6:
                    print("#ERROR: waveform type C parameters insufficient: " + wfp)
                    sys.exit(1)               
                else:
                    I_curr = float( wfp[2])
                    I_floor = float( wfp[3])
                    numOfConsecutiveClk = int( wfp[4])
                    numOfSkippedClk = int( wfp[5])
                    self.AddClkGatingCurr(numOfUnit, I_curr, I_floor, numOfConsecutiveClk, numOfSkippedClk)
            elif wfp[0] == 'D':
                if len(wfp) < 4:
                    print("#ERROR: waveform type D parameters insufficient: " + wfp)
                    sys.exit(1)
                else:
                    self.src_profile_envelope_fileName = wfp[1]
                    self.src_profile_envelope_time_unit_in_sec = float(wfp[2]) 
                    self.src_profile_envelope_waveform_unit = float(wfp[3])
                    self.AddScalingCurr()

            elif wfp[0] == 'E':
                if len(wfp) < 5:
                    print('#ERROR: waveform type E parameters insufficient: ' + wfp)
                    sys.exit(1)
                else:
                    I_floor = float(wfp[2])
                    I_bd_lo = float(wfp[3])
                    I_bd_up = float(wfp[4])
                    if I_bd_lo > I_bd_up:
                        tmp = I_bd_lo
                        I_bd_lo = I_bd_up
                        I_bd_up = tmp 
                    self.AddRandWithinRange(numOfUnit, I_floor, I_bd_lo, I_bd_up)


    def WriteWaveform(self, fileName):
        fout = open(fileName, 'w+')
        fout.write('Time(ns) \t Current(Amp)\n')
        leng_rcd = len( self.currWaveform_list_time_ns)
        for i in range(0, leng_rcd):
            fout.write(str(self.currWaveform_list_time_ns[i]) + '\t' + str(self.currWaveform_list_curr_Amp[i]) + '\n')  
        fout.close()
        print('#INFO: Current waveform output to ' + fileName)

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
print('\n#INFO: Start to run waveform generation at '+ str(voltage) + ' V\n')

waveformInst = CurrWaveform(file_in_para)
waveformInst.CompositeWaveform()
waveformInst.WriteWaveform(file_out_waveform)
waveformInst.WriteWaveform_InTimFormat(file_out_waveform)

if is_print_wf:
    waveformInst.PlotWaveform()

print ("#INFO: Waveform generation exit normally")