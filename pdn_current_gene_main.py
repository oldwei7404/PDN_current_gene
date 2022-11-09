# This script is intended for PDN current profile generation, manipulation 


import os, sys, getopt
import math, cmath
import shutil
import matplotlib.pyplot as plt
import scipy.interpolate

file_dir = ""
file_in_para = ""
file_out_waveform = ""
waveform_params_list = []

###
class CurrWaveform:
    currWaveform_list_time_ns = []  # unit: ns
    currWaveform_list_curr_Amp = []   # unit: Amp
    waveform_params_list = []

    clk_freq = 0
    T_clk = 0
    T_clk_in_ns = 0
    clk_duty_cycle = 0.5
    t_rise_ratio_T = 0.1
    t_fall_ratio_T = 0.1
    nominal_I = 1.0
    I_lkg = 1.0

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
                    cln_str = fin.readline()
                elif cln_str[0] == 'CLK_DutyCycle':
                    self.clk_duty_cycle = float( cln_str[1])
                    print('#INFO: CLK duty cycle:\t' + str(self.clk_duty_cycle))
                    cln_str = fin.readline()
                elif cln_str[0] == 'CLK_T_RISE_as_ratio_of_CLK_Freq':
                    self.t_rise_ratio_T = float( cln_str[1])
                    print('#INFO: CLK rise time (ratio of T):\t' + str(self.t_rise_ratio_T))
                    cln_str = fin.readline()
                elif cln_str[0] == 'CLK_T_FALL_as_ratio_of_CLK_Freq':
                    self.t_fall_ratio_T = float( cln_str[1])
                    print('#INFO: CLK fall time (ratio of T):\t' + str(self.t_rise_ratio_T))
                    cln_str = fin.readline()
                else:  # read in waveform params
                    self.waveform_params_list.append(cln_str_src)
                    wf_cnt = wf_cnt + 1
                    print('INFO: Reading waveform def. ' + str(wf_cnt) + ': ' + self.waveform_params_list[-1])
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

    ### Function: compose the actual waveform based on parameters
    def CompositeWaveform(self):
        for wfp in self.waveform_params_list:
            wfp = wfp.split()
            
            time_wf_ns = float( wfp[1])
            numOfUnit =  int(time_wf_ns/ self.T_clk_in_ns)
        
            if wfp[0] == 'A':
                I_curr = float( wfp[2])
                I_floor = float( wfp[3])
                self.AddConstCLK(numOfUnit, I_curr, I_floor)
            elif wfp[0] == 'B':
                I_start = float( wfp[2])
                I_end = float( wfp[3])
                I_floor = float( wfp[4])          
                self.AddLinearSlopeCurr(numOfUnit, I_start, I_end, I_floor)
            elif wfp[0] == 'C':
                I_curr = float( wfp[2])
                I_floor = float( wfp[3])
                numOfConsecutiveClk = int( wfp[4])
                numOfSkippedClk = int( wfp[5])
                self.AddClkGatingCurr(numOfUnit, I_curr, I_floor, numOfConsecutiveClk, numOfSkippedClk)

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
        print('#INFO: Current waveform TIM format output to ' + fileName)

# Main function 
try:
	opts,args = getopt.getopt(sys.argv[1:],'d:i:o:')
except getopt.GetoptError:
	print('\nUsage: python pdn_current_gene_main.py [-d file directory] [-i input.params] [-o out_curr_profile.tim]')
	sys.exit(2)
if (not opts) and args:
	print('\nUsage: python pdn_current_gene_main.py [-d file directory] [-i input.params] [-o out_curr_profile.tim]')
	sys.exit(2)

for o,a in opts:
    if o=='-d':
        if os.name == 'nt':
            file_dir = a.lstrip(' ').rstrip(' ') + '\\'
        else:
            file_dir = a.lstrip(' ').rstrip(' ') + '/'
    if o=='-i':
        file_in_para = a.lstrip(' ').rstrip(' ')
    if o=='-o':
        file_out_waveform = a.lstrip(' ').rstrip(' ')

# BEGIN read input parameters
file_in_para = file_dir + file_in_para
file_out_waveform = file_dir + file_out_waveform

if os.path.exists(file_in_para):
    print('INFO: input parameters file:    ' + file_in_para)
    print('INFO: output waveform TIM file: ' + file_out_waveform)
else:
    print('ERROR: intput parameter file <' + file_in_para + '> does not exist !')
    sys.exit(1)

# END read input parameters

waveformInst = CurrWaveform(file_in_para)
waveformInst.CompositeWaveform()
waveformInst.WriteWaveform(file_out_waveform)
waveformInst.WriteWaveform_InTimFormat(file_out_waveform)

### interp1d should be enough: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d
### unstructured interpolation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html#scipy.interpolate.griddata



