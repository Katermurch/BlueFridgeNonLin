readout_dict = {
    "RO_LO": 6.6247,
    "RO_LO_pwr": 16,
    "ro_dur": 5000,
}
swap_gate = {
    "swap_freq": 0,
    "swap_amp": 0,
    "swap_time": 0,
}
q1_dict = {
    "qubit_id": "q1",
    "ro_freq": 6.72739,
    "ge_ssm": -0.1144,
    "ef_ssm": -0.2568,
    "ge_time": 82,
    "ef_time": 49,
    "ge_amp": 1,#.5
    "ef_amp":1.48,
    "IQ_angle":110,
    "ro_amp": .15,
    "qubit_thr": [-10000, -600],
}
q2_dict = {
    "qubit_id": "q2",
    "ro_freq": 6.65554,
    "ge_ssm": -0.154,
    "ef_ssm": -0.2962,
    "ge_time": 54,
    "ge_amp": 0.8,
    "ef_amp": 1,
    "IQ_angle": 25,
    "ro_amp": .7,
    "qubit_thr": [-10000, 1900],
}
q3_dict = {
    "qubit_id": "q3",
    "ro_freq": 6.65839,
    "ge_amp": 1,
    "ef_amp": 1,
    "ro_amp": 0.4,
    "RO_LO": 6.4804,
}
q4_dict = {
    "qubit_id": "q4",
    "ro_freq": 6.5113,
    "ge_amp": 1,
    "ef_amp": 1,
    "ro_amp": 0.4,
}
q5_dict = {
    "qubit_id": "q5",
    "ro_freq": 6.44436,
    "ge_amp": 1,
    "ef_amp": 1,
    "ro_amp": 0.4,
}

general_vals_dict = {
    "mixer_offset": 0,
    "mixer_offset_ef": 20,
    "wx_amps": [1,1,1.7,1],  # maximum 1.9
    "coupler_off_value": 0.7,
    "wx_offs": [-0.023, 0, -0.08, -0.096],
    "qubit_bnc": 4.6,
    "TWPA_freq": 4.5,#4.45,
    "TWPA_pwr": -5.4#-5.5,
}
bnc_address = {
    "target_bnc_black": "GPIB0::19::INSTR",
    "big_agilent": "GPIB0::30::INSTR",
    "agilent_function_generator": "GPIB0::30::INSTR",
    "target_bnc_6": "USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR",
    "wx_address": "10.225.208.204",
}
