import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import time
import tewx
SAMPLE_RATE = 1E9 # Gig samples/sec


class Pulse:
    '''
    DESCRIPTION: an object to contain all pulse parameters. 
            If using SSM, self.ssm_bool= True
    
    Parameters for Non-Hermitian qubit    
            if ff is set to any number, the modulation with a cos wave affects the amplitude of the wave
            the t_loop is the period of the sin wave
            the phase_ini is the phase of the sin wave
            (this function made change in gen_pulse) 
    M.A
    PARAMETERS:
        (base): duration, start (time), amplitude, 
        (if ssm_bool): ssm_freq (GHz), phase
    FUNCTIONS:
        make(): create short copy of pulse 
        show(): graph output
        copy(): deep copy, I hope...
        NOTES:
                This does not include difference b/t cos/sin because they can be included in phase
    '''
    def __init__(self,duration, start, amplitude, ssm_freq=None, phase=0, gaussian_bool=False, 
                 phase_ini=0, t_loop=1E-20, ff=None,detunlinear=0.0,detun_NH=0,path_sign=1,detun_NH_phase=0,jmin=0, derived_gaussian_bool = False, ramp_up_bool = False, ramp_down_bool = False):
        self.duration = int(duration)
        self.start = int(start)
        self.amplitude = amplitude
        self.phase_ini = phase_ini
        self.t_loop = t_loop
        self.ff=ff
        self.path_sign=path_sign
        self.detun_NH = detun_NH
        self.detunlinear=detunlinear
        self.detun_NH_phase = detun_NH_phase
        self.jmin=jmin
        if ssm_freq is not None:
            self.ssm_bool = True  
            self.ssm_freq = ssm_freq
            self.phase = phase

        else:
            self.ssm_bool = False
        #self.waveform = self.make() ## make is currently not working.

        # FUTURE FEATURES:
        self.gaussian_bool = gaussian_bool
        self.derived_gaussian_bool = derived_gaussian_bool
        self.ramp_up_bool = ramp_up_bool
        self.ramp_down_bool = ramp_down_bool

    
    def make(self):
        new_array = np.zeros(self.duration)
        if self.ssm_bool:
            gen_pulse( dest_wave = new_array, pulse=self)
        else:
            gen_pulse(new_array, pulse=self)

        return new_array

    def show(self):
        plt.plot(np.arange(self.start,self.duration), self.waveform)
        plt.show()

    def copy(self):
        # there must be a better/ more general way to do this.
        if self.ssm_bool:
            return Pulse( self.duration, self.start, self.amplitude, self.ssm_freq, self.phase, self.gaussian_bool, 
                         self.phase_ini, self.t_loop, self.ff,self.detunlinear,self.detun_NH,self.path_sign,self.detun_NH_phase,self.jmin,self.derived_gaussian_bool, self.ramp_up_bool, self.ramp_down_bool)
        else:
            return Pulse(self.duration, self.start, self.amplitude, gaussian_bool=self.gaussian_bool)

    def toString(self):
        outString = "Pulse of {0} [amp] from {1}+{2}".format(self.amplitude, self.start, self.duration)
        if self.ssm_bool:
            outString += " SSM @ {0} MHz with phase={1}".format(self.ssm_freq*1000, self.phase)
        return outString
#END pulse class        
    
class Sequence:
    '''
    DESCRIPTION: an object to contain all sequence parameters. 
    PARAMETERS:
        (base): duration, start (time), amplitude, 
        (if ssm_bool): ssm_freq (GHz), phase
    FUNCTIONS:
        
        NOTES:
                This does not include difference b/t cos/sin because they can be included in phase
    '''
    def __init__(self,sequence_length, num_steps,mixer_orthogonality=90):
        
        self.sequence_length = int(sequence_length)
        self.num_steps = int(num_steps)
        self.mixer_orthogonality = mixer_orthogonality

        self.channel_list = self._initialize_channels()
        # all_data is list of [ch1, ch2, ch3, ch4]
        #   e.g. ch1 =  [waveform, m1, m2]
        #       e.g. waveform contains sweep: m1.shape = [num_steps, samples_per_step]
        
    def insert_waveform(self, channel_num, pulse, step_index):
        full_seq_waveform = self.all_data[channel_num-1][0] # 0 for waveform
        current_step = full_seq_waveform[step_index]
        gen_pulse(current_step, pulse) ##in-place insertion
        #self.all_data[channel_num-1][0] += current_step
        
    def insert_marker(channel_num, marker_num, pulse):
        full_seq_marker = self.all_data[channel_num-1][marker_num-1]
        current_step = full_seq_marker[step_index]
        gen_pulse(current_step, pulse)
        self.all_data[channel_num-1][marker_num-1] += current_step

    def insert_bothChannels(self,primary_channel_num, pulse,step_index):
        ## adds pulse to both channels, offset by mixer_ortho.
        ch_num = primary_channel_num
        self.insert_waveform(ch_num,pulse,step_index)
        copy = pulse.copy()
        copy.phase += self.mixer_orthogonality
        self.insert_waveform(ch_num,copy,step_index)
        
    def convert_to_tabor_format(self, channel_num):
        # each of the following is a [num_steps x samples_per_step] matrix
        waveform = self.all_data[channel_num][0]
        mark1 = self.all_data[channel_num][1]
        mark2 = self.all_data[channel_num][2]
        binarized = int(2**12*waveform) + int(2**14 *mark1) + int(2**15 *mark2)       
        return binarized
    
    
    def add_gate(self, source_1, source_2=None,destination_tuple=(1,1)): #input channel numbers #channel1 marker1 is default use dest_tuple=(3,2) for ch3/4 mkr2
        # each of the following is a [num_steps x samples_per_step] matrix
        channel1_channel = self.channel_list[source_1-1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
        both_ch1_ch2 = channel1_channel**2
        if source_2:
            channel2_channel = self.channel_list[source_2-1][0] # dim 0: channel 1; dim 1: [ch,m1,m2]
            both_ch1_ch2 += channel2_channel**2
        qubit_gate = create_gate(both_ch1_ch2)
        self.channel_list[destination_tuple[0]-1][destination_tuple[1]] = qubit_gate     

    def add_sweep(self, channel, marker=0, sweep_name='none', start=0, stop=0, initial_pulse=Pulse(amplitude=0, duration=0, start=0)):
        '''
        DESCRIPTION:  A thin wrapper to add pulses to the correct channel 
        INPUT: channel = {1-4} (converts to indices); marker=0(channel),1,2; arguments for gen_sweep 
        OUTPUT:
        '''

        ## error checking
        if start==stop and sweep_name != 'none':
            raise Warning("Start and sweep are the same; did you mean that?")
        if channel not in [1,2,3,4]:
            raise IOError("Invalid channel number: "+str(channel))
        if marker not in [0,1,2]:
            raise IOError("Invalid marker number: "+str(marker))
       
        ## send the input to _gen_sweep 
        dest_wave = self.channel_list[channel-1][marker]
        self._gen_sweep(sweep_name, start=start, stop=stop, dest_wave=dest_wave, initial_pulse=initial_pulse )
    #END add_sweep


    def _initialize_sequence_matrix(self):
        '''
        DESCRIPTION: prepare an empty matrix of size [time_steps, sequence_steps] for each channel
        INPUT:
        OUTPUT:
        '''
        num_steps = self.num_steps
        file_length = self.sequence_length

        channel= np.zeros((num_steps,file_length))
        mark1  = np.zeros((num_steps,file_length))
        mark2  = np.zeros((num_steps,file_length))
        
        return channel, mark1, mark2
    ##END _intiialize_sequence_matrix

    
    def _initialize_channels(self):
        '''
        DESCRIPTION: prepare the channels and markers
        INPUT:
        OUTPUT:
        '''
        num_channels = 4
        
        channel_array = [ [None,None,None] for i in range(num_channels) ] 
        for ch_index in range(len(channel_array)):
            wave, mark1, mark2 = self._initialize_sequence_matrix()
            channel_array[ch_index] = [wave, mark1, mark2]
        
        # The WX2184C channel (ch) 1 and 2 share markers; likewise ch 3 and 4
        #   so we will copy ch 1 to 2 and 3 to 4
        mark1_index, mark2_index = 1,2
        ## set ch2 markers equal to ch1 markers
        channel_array[1][mark1_index] = channel_array[0][mark1_index] # ch1/2 m1
        channel_array[1][mark2_index] = channel_array[0][mark2_index] # ch1/2 m2
        ## set ch4 markers equal to ch3 markers
        channel_array[3][mark1_index] = channel_array[2][mark1_index] # ch3/4 mmark1_index
        channel_array[3][mark2_index] = channel_array[2][mark2_index] # ch3/4 mmark2_index
        return channel_array
    ##END _intitialize_channels

    
    def _gen_sweep(self, sweep_name, start, stop, dest_wave, initial_pulse=Pulse(amplitude=0, duration=0,start=0)):
        '''
        DESCRIPTION:  sweeps 'none', 'amplitude', 'width', 'start' or 'phase' 
                'none' sets same pulse for all steps in sequence
        INPUT:       sets range for initial + [start,stop)
                updates parameters to initial_pulse
        OUTPUT:     writes to dest_wave 
        '''
#        phase_ini = self.phase_ini
#        t_loop=self.t_loop
#        detun=Pulse.detun # this is a detuning based Cos function for encircling EP 
#        detunlin = self.detunlin # this is a fixed detuning
#        path_sign=self.path_sign
#        detun_NH_phase=self.detun_NH_phase
        
        
        ## Check input
        if len(dest_wave) < self.num_steps:
            raise IOError("dest_wave is too short ({0})".format(len(dest_wave)))
        
        updated_pulse = initial_pulse.copy()
        if sweep_name == 'none':
            for step_index in range(self.num_steps):
                gen_pulse( dest_wave = dest_wave[step_index], pulse=updated_pulse)

        elif sweep_name == 'start':

            for step_index, param_val in enumerate(np.linspace(start,stop,self.num_steps)):
                updated_pulse.start = initial_pulse.start  + int(param_val) #- initial_pulse.duration
                gen_pulse( dest_wave = dest_wave[step_index], pulse=updated_pulse)
                
        elif sweep_name == 'amplitude':
            for step_index, param_val in enumerate(np.linspace(start,stop,self.num_steps)):
                updated_pulse.amplitude= initial_pulse.amplitude+ param_val
                gen_pulse( dest_wave = dest_wave[step_index], pulse=updated_pulse)

        elif sweep_name == 'width':
            for step_index, param_val in enumerate(np.linspace(start,stop,self.num_steps)):
                updated_pulse.duration= initial_pulse.duration + int(param_val)
                gen_pulse( dest_wave = dest_wave[step_index], pulse=updated_pulse)

        elif sweep_name == 'phase':
            if not initial_pulse.ssm_bool: raise ValueError("Sweeping phase w/o SSM")
            for step_index, param_val in enumerate(np.linspace(start,stop,self.num_steps)):
                updated_pulse.phase= initial_pulse.phase+ param_val
                gen_pulse( dest_wave = dest_wave[step_index], pulse=updated_pulse)
                
        elif sweep_name == 'phase_linear_detun':
            if not initial_pulse.ssm_bool: raise ValueError("Sweeping phase w/o SSM")
            for step_index, param_val in enumerate(np.linspace(start,stop,self.num_steps)):
                updated_pulse.phase= initial_pulse.phase - param_val*initial_pulse.detunlinear*360
                gen_pulse( dest_wave = dest_wave[step_index], pulse=updated_pulse)
                
        elif sweep_name == 'detuned_phase':
            if not initial_pulse.ssm_bool: raise ValueError("Sweeping phase w/o SSM")
            for step_index, param_val in enumerate(np.linspace(start,stop,self.num_steps)):
#                - 1/(1/t_loop)*detun_NH*1E9/SAMPLE_RATE*np.cos(2*np.pi*(times)/t_loop+phase_ini)
                updated_pulse.phase= initial_pulse.phase+ (initial_pulse.detun_NH_phase*initial_pulse.t_loop*np.cos(2*np.pi*param_val/initial_pulse.t_loop + initial_pulse.phase_ini)-initial_pulse.detun_NH_phase*initial_pulse.t_loop)*180/np.pi
                gen_pulse( dest_wave = dest_wave[step_index], pulse=updated_pulse)
        elif sweep_name == 'Jmin_var':
            initial_pulse.jmin=0
            for step_index, param_val in enumerate(np.linspace(start,stop,self.num_steps)):
#                - 1/(1/t_loop)*detun_NH*1E9/SAMPLE_RATE*np.cos(2*np.pi*(times)/t_loop+phase_ini)
                updated_pulse.jmin= initial_pulse.jmin+ param_val
                gen_pulse( dest_wave = dest_wave[step_index], pulse=updated_pulse)

        else:
            raise ValueError("Bad sweep parameter: "+sweep_name)
    #END gen_sweep 


    def write_sequence(self, base_name='foo', file_path=os.getcwd(), use_range_01=False,num_offset=0, write_binary=False, pr = 1): 
        ''' 
        DESCRIPTION: writes a single channel INPUT:
            (optional) mark1/2: numpy arrays with marker data (0,1)
        OUTPUT:
        TODO:
        '''
        
        if not file_path.endswith("\\"): file_path+= "\\"
        print("writing to {}".format(file_path))

        for ch_index, (channel, mark1, mark2) in enumerate(self.channel_list):
            ch_name = "ch" + str(ch_index+1)
            if pr==1:
                print("writing "+ch_name)
            
            for step_index in range(self.num_steps):
                if write_binary:
                    file_name = file_path+base_name+"_"+ch_name+"_{:04d}.npy".format(step_index+num_offset)
                else:
                    file_name = file_path+base_name+"_"+ch_name+"_{:04d}.csv".format(step_index+num_offset)

                if use_range_01: # write floats between 0 and 1
                    with_index = zip(range(len(channel[step_index])), channel[step_index] )
                    np.savetxt(file_name, with_index, fmt='%d, %f')
                    continue

                # convert to binary
                # 15 bits = 12 bits of information, 2**13=sign bit, 2**14=mark1, 2**15=mark2
                else:
                    new_mark1 = [ int(i*2**14) for i in mark1[step_index] ]
                    new_mark2 = [ int(i*2**15) for i in mark2[step_index] ]
                    if write_binary:
                        binary_file = np.array([ round(2**12 *val + 8191.5) for val in channel[step_index] ])
                        binary_file = binary_file.clip(0, 2**14
                                                       -1)
                    else:
                        binary_file = np.array([ int(2**12 *val) for val in channel[step_index] ])
                    binary_file += new_mark1
                    binary_file += new_mark2
                    
                    #pd.DataFrame(binary_file).to_csv(file_name,float_format='%d', header=False)
                    if write_binary:
                        binary_file = binary_file.astype('uint16')
                        np.save(file_name, binary_file)
                    else:
                        with_index = zip(range(len(binary_file)), binary_file)
                        np.savetxt(file_name, list(with_index), fmt='%d, %d')
            ##end for loop through channels
    #END write_sequence
    
    def load_sequence(self, instr_addr, base_name='foo', file_path=os.getcwd(), num_offset=0,amp12=1.5,amp34=1.5):
        ''' 
        DESCRIPTION: loads multi channel INPUT:
        OUTPUT:
        TODO:
        '''
        file_length = self.sequence_length
        num_steps = self.num_steps
        if not file_path.endswith("\\"): file_path+= "\\"
        print("loading {}".format(file_path))
 
        # Reading the wave data from .npy file
        waveforms = [[ None, ] * num_steps for _ in self.channel_list]
        for ch_index, _ in enumerate(self.channel_list):
            ch_name = "ch" + str(ch_index+1)
            print("loading "+ch_name)
            for step_index in range(num_steps):
                file_name = file_path+base_name+"_"+ch_name+"_{:04d}.npy".format(step_index+num_offset)
                waveforms[ch_index][step_index] = np.load(file_name)
       
        # Initializing the instrument
        inst = tewx.TEWXAwg(instr_addr, paranoia_level=1)
        inst.send_cmd('*CLS') # Clear errors
        inst.send_cmd('*RST') # Reset the device #need to add several commands to set up device to use markers and other configurations
        inst.send_cmd(':OUTP:ALL 0')
        seg_quantum = inst.get_dev_property('seg_quantum', 16)
        
        # Setting up the markers
        
        
        # Downloading the wave data
        seg_len = np.ones(num_steps, dtype=np.uint32) * file_length
        pseudo_seg_len = num_steps * file_length + (num_steps - 1) * seg_quantum
        wav_dat = np.zeros(2 * pseudo_seg_len, 'uint16')   
  
        for ch_index, _ in enumerate(self.channel_list):
            if ch_index % 2:
                continue
            offs = 0
            for step_index in range(self.num_steps):
                wav1 = waveforms[ch_index][step_index]
                wav2 = waveforms[ch_index + 1][step_index]
                offs = inst.make_combined_wave(wav1, wav2, wav_dat, dest_array_offset=offs, add_idle_pts=(0!=offs))
            
            # select channel:
            inst.send_cmd(':INST:SEL {0}'.format(ch_index+1))
            
            inst.send_cmd('MARK:SEL 1')
            inst.send_cmd('MARK:SOUR USER')
            inst.send_cmd('MARK:STAT ON')
            inst.send_cmd('MARK:SEL 2')
            inst.send_cmd('MARK:SOUR USER')
            inst.send_cmd('MARK:STAT ON')
            # select user-mode (arbitrary-wave):
            inst.send_cmd(':FUNC:MODE FIX')
            # delete all segments (just to be sure):
            inst.send_cmd(':TRAC:DEL:ALL')
            inst.send_cmd('SEQ:DEL:ALL')
            # set combined wave-downloading-mode:
            inst.send_cmd(':TRAC:MODE COMB')
            # define the pseudo segment:
            inst.send_cmd(':TRAC:DEF 1,{0}'.format(np.uint32(pseudo_seg_len)))
            # select segment 1:
            inst.send_cmd(':TRAC:SEL 1')
            # download binary data:
            
            inst.send_binary_data(':TRAC:DATA', wav_dat)
            
            # ---------------------------------------------------------------------
            # Write the *appropriate* segment-table
            # (array of 'uint32' values holding the segments lengths)
            # ---------------------------------------------------------------------
            inst.send_binary_data(':SEGM:DATA', seg_len)
            # Setting up sequence mode
            for step in range(1, num_steps + 1):
                inst.send_cmd(':SEQ:DEF {},{},1,0'.format(step, step))
            inst.send_cmd(':FUNC:MODE SEQ')
            inst.send_cmd(':SEQ:ADV STEP')
            # Setting up the triggers
            inst.send_cmd(':TRIG:SOUR EVEN')
            inst.send_cmd(':TRIG:COUN 1')
            # Turn channels on:
            inst.send_cmd(':INIT:CONT 0')
            
        # Setting up amplitudes and offsets
        amp = [amp12, amp12, amp34,amp34]
        offset = [0., 0., 0.057, -0.061]
        for ch_index, _ in enumerate(self.channel_list):
            inst.send_cmd(':INST:SEL {0}'.format(ch_index+1))
            inst.send_cmd(':VOLT {}'.format(amp[ch_index]))
            inst.send_cmd(':VOLT:OFFS {}'.format(offset[ch_index]))
            

        
        inst.send_cmd(':INST:COUP:STAT ON')
        inst.send_cmd(':OUTP:ALL 1')
        

        

        
        
        # query system error
        syst_err = inst.send_query(':SYST:ERR?')
        print(syst_err)
        inst.close()
        
    def load_sequence_CSV(self, instr_addr, base_name='foo', file_path=os.getcwd(), num_offset=0):
        ''' 
        DESCRIPTION: loads multi channel INPUT:
        OUTPUT:
        TODO:
        '''
        file_length = self.sequence_length
        num_steps = self.num_steps
        if not file_path.endswith("\\"): file_path+= "\\"
        print("loading {}".format(file_path))
 
        # Reading the wave data from .npy file
        waveforms = [[ None, ] * num_steps for _ in self.channel_list]
        for ch_index, _ in enumerate(self.channel_list):
            ch_name = "ch" + str(ch_index+1)
            print("loading "+ch_name)
            for step_index in range(num_steps):
                file_name = file_path+base_name+"_"+ch_name+"_{:d}.csv".format(step_index+num_offset)
                waveforms[ch_index][step_index] = np.loadtxt(file_name,dtype=int,delimiter=', ',usecols=(1,))
                
       
        # Initializing the instrument
        inst = tewx.TEWXAwg(instr_addr, paranoia_level=1)
        inst.send_cmd('*CLS') # Clear errors
        inst.send_cmd('*RST') # Reset the device #need to add several commands to set up device to use markers and other configurations
        inst.send_cmd(':OUTP:ALL 0')
        seg_quantum = inst.get_dev_property('seg_quantum', 16)
        
        # Setting up the markers
        
        
        # Downloading the wave data
        seg_len = np.ones(num_steps, dtype=np.uint32) * file_length
        pseudo_seg_len = num_steps * file_length + (num_steps - 1) * seg_quantum
        wav_dat = np.zeros(2 * pseudo_seg_len, 'uint16')   
  
        for ch_index, _ in enumerate(self.channel_list):
            if ch_index % 2:
                continue
            offs = 0
            for step_index in range(self.num_steps):
                wav1 = waveforms[ch_index][step_index]
                wav2 = waveforms[ch_index + 1][step_index]
                offs = inst.make_combined_wave(wav1, wav2, wav_dat, dest_array_offset=offs, add_idle_pts=(0!=offs))
            
            # select channel:
            inst.send_cmd(':INST:SEL {0}'.format(ch_index+1))
            
            inst.send_cmd('MARK:SEL 1')
            inst.send_cmd('MARK:SOUR USER')
            inst.send_cmd('MARK:STAT ON')
            inst.send_cmd('MARK:SEL 2')
            inst.send_cmd('MARK:SOUR USER')
            inst.send_cmd('MARK:STAT ON')
            # select user-mode (arbitrary-wave):
            inst.send_cmd(':FUNC:MODE FIX')
            # delete all segments (just to be sure):
            inst.send_cmd(':TRAC:DEL:ALL')
            inst.send_cmd('SEQ:DEL:ALL')
            # set combined wave-downloading-mode:
            inst.send_cmd(':TRAC:MODE COMB')
            # define the pseudo segment:
            inst.send_cmd(':TRAC:DEF 1,{0}'.format(np.uint32(pseudo_seg_len)))
            # select segment 1:
            inst.send_cmd(':TRAC:SEL 1')
            # download binary data:
            inst.send_binary_data(':TRAC:DATA', wav_dat)
            
            # ---------------------------------------------------------------------
            # Write the *appropriate* segment-table
            # (array of 'uint32' values holding the segments lengths)
            # ---------------------------------------------------------------------
            inst.send_binary_data(':SEGM:DATA', seg_len)
            # Setting up sequence mode
            for step in range(1, num_steps + 1):
                inst.send_cmd(':SEQ:DEF {},{},1,0'.format(step, step))
            inst.send_cmd(':FUNC:MODE SEQ')
            inst.send_cmd(':SEQ:ADV STEP')
            # Setting up the triggers
            inst.send_cmd(':TRIG:SOUR EVEN')
            inst.send_cmd(':TRIG:COUN 1')
            # Turn channels on:
            inst.send_cmd(':INIT:CONT 0')
            
        # Setting up amplitudes and offsets
        amp = [1, 1.02, 0.5, .5]
        offset = [-0.04, -0.053, -0.013, -0.002]
        for ch_index, _ in enumerate(self.channel_list):
            inst.send_cmd(':INST:SEL {0}'.format(ch_index+1))
            inst.send_cmd(':VOLT {}'.format(amp[ch_index]))
            inst.send_cmd(':VOLT:OFFS {}'.format(offset[ch_index]))
            

        
        inst.send_cmd(':INST:COUP:STAT ON')
        inst.send_cmd(':OUTP:ALL 1')
        

        

        
        
        # query system error
        syst_err = inst.send_query(':SYST:ERR?')
        print(syst_err)
        inst.close()     
                
        
    
    def convert_to_tabor_format(self):
        '''
        DESCRIPTION: converts the sequence structure to a loadable format
        INPUT: populated Sequence
        OUTPUT: array of size [4, num_steps, sequence_length] in binarized form with markers
        '''
        tabor_format = np.zeros((4,self.num_steps, self.sequence_length))
        for ch_index, (channel, mark1, mark2) in enumerate(self.channel_list):
            ## loop through ch 1-4
            single_channel_binary  = np.array(2**12 *channel,dtype='int')
            single_channel_binary += np.array(2**14 *mark1,dtype='int')
            single_channel_binary += np.array(2**15 *mark2,dtype='int')
            tabor_format[ch_index] = single_channel_binary
        ##END loop through channels
        return tabor_format
    ##END convert_to_tabor_format
##END Sequence 


def gen_pulse(dest_wave, pulse):
    ## consider renaming to insert_pulse()
    
    '''
    DESCRIPTION: generates pulse on one wave
            Note, this does not add constant pulese through all steps; that's handled by add_sweep('none')
    INPUT: 
        Pulse object (contains start,duration, etc)
    OUTPUT:
        in-place adjustment to dest_wave
    NOTES:
    TODO:
    '''
    ## Decompose pulse object
    start = pulse.start
    dur = pulse.duration
    amp = pulse.amplitude
    phase_ini = pulse.phase_ini*np.pi/180
    t_loop=pulse.t_loop
    ff=pulse.ff
#    detun=pulse.detun
    path_sign=pulse.path_sign
    detunlinear=pulse.detunlinear
    detun_NH=pulse.detun_NH
    detun_NH_phase=pulse.detun_NH_phase
    bmin=1+pulse.jmin

    if dur <0:
        dur = abs(dur)
        start -= dur
    if ff==None:
          ##  Create output
          if pulse.ssm_bool:
              ssm_freq = pulse.ssm_freq
              phase = pulse.phase
          
              # start times depend on start and duration because curves should pick up same absolute phase( may be shifted for cos/sin/etc), ie two pulses placed front to back should continue overall
              times = np.arange(start,start+dur)
              ang_freq = 2*np.pi*(ssm_freq*1E9)/SAMPLE_RATE # convert to units of SAMPLE_RATE
              phase_rad = phase/180.0*np.pi #- 2*np.pi*detunlinear*(start) # - 2*np.pi*detun_NH_phase*np.sin(2*np.pi*(times)/t_loop+phase_ini)#
              addition = amp*np.sin(ang_freq*times + phase_rad)
          else: 
              addition = amp* np.ones(dur)    
              
          if pulse.gaussian_bool:
              argument = -(times-start-dur/2)**2
              argument /= 2*(dur*0.2)**2 # 0.847 gives Gauss(start +dur/2) = 0.5
              gauss_envelope = np.exp(argument);
              addition *= gauss_envelope  
      
          if pulse.derived_gaussian_bool:
              argument = -(times-start-dur/2)**2
              argument /= 2*(dur*0.2)**2 # 0.847 gives Gauss(start +dur/2) = 0.5
              gauss_envelope = np.exp(argument);
              
              derivative = -(times-start-dur/2)
              derivative /= (dur*0.2)**2
              derived_gaussian = derivative*gauss_envelope
              addition *= derived_gaussian
              
          if pulse.ramp_up_bool:
              rampUpFunc = 1/dur*(times-start)
              addition *= rampUpFunc
              
          if pulse.ramp_down_bool:
              rampDownFunc = 1 - 1/dur*(times-start)
              addition *= rampDownFunc
              
    elif ff==1:
          if pulse.ssm_bool:
              ssm_freq = pulse.ssm_freq    
              phase = pulse.phase
          
              # start times depend on start and duration because curves should pick up same absolute phase( may be shifted for cos/sin/etc), ie two pulses placed front to back should continue overall
              times = np.arange(start,start+dur)
              ang_freq = 2*np.pi*(ssm_freq*1E9)/SAMPLE_RATE # convert to units of SAMPLE_RATE
              phase_rad = phase/180.0*np.pi
              ampfunc = 1/2*((2-bmin)*np.cos(2*np.pi*(times-start)/t_loop+phase_ini)+bmin)
#              freqfunc= np.sin()
#              addition = ampfunc*amp*np.sin(ang_freq*times + phase_rad)
              omega=2*np.pi/t_loop
#              freqfunc= 2*np.pi*(detun*1E-3*1E9)/SAMPLE_RATE
              addition = ampfunc*amp*np.sin((t_loop*detun_NH*1E9/SAMPLE_RATE*np.cos(2*np.pi*(times-start)/t_loop+phase_ini)-t_loop*detun_NH*1E9/SAMPLE_RATE)+ang_freq*times + phase_rad) #frequency modulation through phase                                             
             
#          if pulse.ssm_bool:
#              ssm_freq = pulse.ssm_freq    
#              phase = pulse.phase
#              
#              #Begin by ramping up the J to its initial value
#              pre_times = np.arange(start, start+j_ramp_time)
#              ang_freq = 2*np.pi*(ssm_freq*1E9)/SAMPLE_RATE # convert to units of SAMPLE_RATE
#              phase_rad = phase/180.0*np.pi
#              rampUpFunc = amp/j_ramp_time*(pre_times-pre_times[0])
#              addition = rampUpFunc*np.sin(1/(1/t_loop)*detun_NH*1E9/SAMPLE_RATE*np.cos(2*np.pi*(pre_times-start)/t_loop+phase_ini)-1/(1/t_loop)*detun_NH*1E9/SAMPLE_RATE+ang_freq*pre_times + phase_rad)
#
#              # start times depend on start and duration because curves should pick up same absolute phase( may be shifted for cos/sin/etc), ie two pulses placed front to back should continue overall
#              times = np.arange(start+j_ramp_time, start+dur-j_ramp_time)
#              ang_freq = 2*np.pi*(ssm_freq*1E9)/SAMPLE_RATE # convert to units of SAMPLE_RATE
#              phase_rad = phase/180.0*np.pi
#              ampfunc = 1/2*((2-bmin)*np.cos(2*np.pi*(times-times[0])/(dur)+phase_ini)+bmin)
##              freqfunc= np.sin()
##              addition = ampfunc*amp*np.sin(ang_freq*times + phase_rad)
#              omega=2*np.pi/t_loop
##              freqfunc= 2*np.pi*(detun*1E-3*1E9)/SAMPLE_RATE
#              addition += ampfunc*amp*np.sin(1/(1/t_loop)*detun_NH*1E9/SAMPLE_RATE*np.cos(2*np.pi*(times-start)/t_loop+phase_ini)-1/(1/t_loop)*detun_NH*1E9/SAMPLE_RATE+ang_freq*times + phase_rad)
#                  
#              #End by ramping J down to 0
#              post_times = np.arange(start+dur-j_ramp_time, start+dur)
#              ang_freq = 2*np.pi*(ssm_freq*1E9)/SAMPLE_RATE # convert to units of SAMPLE_RATE
#              phase_rad = phase/180.0*np.pi
#              rampDownFunc = amp - amp/j_ramp_time*(post_times-post_times[0])
#              addition += rampDownFunc*np.sin(1/(1/t_loop)*detun_NH*1E9/SAMPLE_RATE*np.cos(2*np.pi*(post_times-start)/t_loop+phase_ini)-1/(1/t_loop)*detun_NH*1E9/SAMPLE_RATE+ang_freq*post_times + phase_rad)
              
          else: 
              addition = amp* np.ones(dur)    
              
          if pulse.gaussian_bool:
              argument = -(times-start-dur/2)**2
              argument /= 2*(dur*0.2)**2 # 0.847 gives Gauss(start +dur/2) = 0.5
              gauss_envelope = np.exp(argument);
              addition *= gauss_envelope  
      
          if pulse.derived_gaussian_bool:
              argument = -(times-start-dur/2)**2
              argument /= 2*(dur*0.2)**2 # 0.847 gives Gauss(start +dur/2) = 0.5
              gauss_envelope = np.exp(argument);
              
              derivative = -(times-start-dur/2)
              derivative /= (dur*0.2)**2
              derived_gaussian = derivative*gauss_envelope
              addition *= derived_gaussian
              
          if pulse.ramp_up_bool:
              rampUpFunc = amp/dur*(times-start)
              addition *= rampUpFunc
              
          if pulse.ramp_down_bool:
              rampDownFunc = amp - amp/dur*(times-start)
              addition *= rampDownFunc
          
          
          
    try: 
        dest_wave[start:start+dur] += addition
    except ValueError:
        print( "Over-extended pulse (ignored):\n {0}".format(pulse.toString()))
#END gen_pulse


def some_Fun():
    '''
    DESCRIPTION: 
    INPUT:
    OUTPUT:
    TODO:
    '''
    pass

'''
def rabi_seq():
    file_length= 8000 # becomes global FILELENGTH
    num_steps = 101# becomes global NUMSTEPS
    #WAVE,MARK1, MARK2 = initialize(file_length, num_steps)
    ALL_CHANNELS = initialize_wx(file_length, num_steps)
    # ALL_CHANNELS is 4-array for Ch 1,2. Each elem is a tuple of (channel, M1,M2)
    #       each tuple elem. is a seq_len X num_samples matrix.

    ## channels
    p = Pulse(start=5795, duration=0,  amplitude=0.5, ssm_freq=0.200, phase=0)
    add_sweep(ALL_CHANNELS, channel=1, sweep_name='width', start=0 , stop= 200, initial_pulse=p )
    p.phase = 98.0
    add_sweep(ALL_CHANNELS, channel=2, sweep_name='width', start=0 , stop= 200, initial_pulse=p )
   
    readout = Pulse(start=6000,duration=1000,amplitude=1)
    add_sweep(ALL_CHANNELS, channel=3, sweep_name='none',initial_pulse=readout)

    ## markers
    gate = Pulse(start=5790, duration=10, amplitude=1)
    add_sweep(ALL_CHANNELS, channel=1, marker=1, sweep_name='width', start=0, stop=220, initial_pulse=gate)
    trigger = Pulse(start=2000, duration=1000, amplitude=1)
    add_sweep(ALL_CHANNELS, channel=3, marker=1, sweep_name='none', initial_pulse=trigger)
    
    return ALL_CHANNELS
    ## send to ARB
    dir_name = r"C:\Arb Sequences\EUR_sequences\mixerOrthogonality_98deg\piTime_23ns\tmp"
    write_sequence(ALL_CHANNELS, file_path=dir_name, use_range_01=False)
##END rabi_seq
'''


def create_gate(seq_matrix,width=5):
    '''
    DESCRIPTION: for all times (ts) and for all steps (ss): if any amplitude exists, extend in time by width
    INPUT: seq_matrix of size [sequence_steps, samples]
    OUTPUT: binary mask with same size as input
    ## KNOWN BUG: if pulse_end + width > num_samples will create error.
    '''
    mask_ss, mask_ts= np.where(seq_matrix != 0)
    gate = seq_matrix.copy()
    gate[ (mask_ss, mask_ts)] = 1
    gate[ (mask_ss, mask_ts-width)] = 1
    gate[ (mask_ss, mask_ts+width)] = 1
   
    return gate
##END create_gate
