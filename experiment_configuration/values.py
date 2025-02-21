readout_dict = {
    "RO_LO": 6.6247,
    "RO_LO_pwr": 16,
    "ro_dur": 4000,
}
q1_dict = {
    "qubit_id": "q1",
    "ro_freq": 6.727420,
    "ge_ssm": -0.114,
    "ef_ssm": -0.186,
    "ge_time":  66.625,
    "ef_time": 77,
    "ge_amp": .5,
    "ef_amp":1,
    "IQ_angle": 205,
    "ro_amp": 0.6,
    "qubit_thr": [-10000,-600],
}
q2_dict = {
    "qubit_id": "q2",
    "ro_freq": 6.6556,
    "ge_ssm": -0.154,
    "ef_ssm": -0.224,
    "ge_time": 40.15374412400174,
    "ge_amp": 1,
    "ef_amp": 1,
    "IQ_angle": 85,
    "ro_amp": 0.4,
    "qubit_thr": [-10000,1900],
}
q3_dict = {
    "qubit_id": "q3",
    "ro_freq": 6.65839,
    "ge_amp": 1,
    "ef_amp": 1,
    "ro_amp": 0.4,
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
    "wx_amps": [1.0, 1, 0.5, 1],
    "coupler_off_value": 0.7,
    "wx_offs": [0.7, 0, 0, 0],
    "qubit_bnc": 4.6,
    "TWPA_freq": 4.5,
    "TWPA_pwr": -5.4
}
bnc_address = {
    "target_bnc_black": "GPIB0::19::INSTR",
    "big_agilent": "GPIB0::30::INSTR",
    "agilent_function_generator": 'GPIB0::30::INSTR',
    "target_bnc_6":"USB0::0x03EB::0xAFFF::411-433500000-0753::INSTR",
    "wx_address": '128.252.134.31'

}
