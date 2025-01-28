class Qubit:
    def __init__(
        self,
        qubit_id,
        ro_freq=None,
        ge_time=None,
        ef_time=None,
        ge_ssm=None,
        ef_ssm=None,
        ge_amp=None,
    ):
        """
        Initialize a Qubit instance with specific properties.

        Args:
            qubit_id (str): The identifier for the qubit (e.g., "q1", "q2").
            ro_freq (float): Readout frequency for the qubit.
            ge_time (int): Gate time for the ground-to-excited (ge) state transition.
            ef_time (int): Gate time for the excited-to-first (ef) state transition.
            ge_ssm (float): Sideband modulation for ge transition.
            ef_ssm (float): Sideband modulation for ef transition.
            ge_amp (float): Amplitude for ge transition.
        """
        self.qubit_id = qubit_id
        self.ro_freq = ro_freq
        self.ge_time = ge_time
        self.ef_time = ef_time
        self.ge_ssm = ge_ssm
        self.ef_ssm = ef_ssm
        self.ge_amp = ge_amp

    def __repr__(self):
        return (
            f"Qubit({self.qubit_id}, ro_freq={self.ro_freq}, ge_time={self.ge_time}, "
            f"ef_time={self.ef_time}, ge_ssm={self.ge_ssm}, ef_ssm={self.ef_ssm}, ge_amp={self.ge_amp})"
        )

class Readout:
    def __init__(
        self,
        readout_amp_1=None,
        readout_amp_2=None,
        RO_LO=None,
        RO_LO_pwr=None,
        ROq3=None,
        ro_pulse_duration=None,
        ro_dur=None,
    ):
        """
        Initialize a Qubit instance with specific properties.

        Args:
            qubit_id (str): The identifier for the qubit (e.g., "q1", "q2").
            ro_freq (float): Readout frequency for the qubit.
            ge_time (int): Gate time for the ground-to-excited (ge) state transition.
            ef_time (int): Gate time for the excited-to-first (ef) state transition.
            ge_ssm (float): Sideband modulation for ge transition.
            ef_ssm (float): Sideband modulation for ef transition.
        """
        self.readout_amp_1 = readout_amp_1
        self.readout_amp_2 = readout_amp_2
        self.RO_LO = RO_LO
        self.RO_LO_pwr = RO_LO_pwr
        self.ROq3 = ROq3
        self.ro_pulse_duration = ro_pulse_duration
        self.ro_dur = ro_dur

    def __repr__(self):
        return (
            f"Readout(readout_amp_1={self.readout_amp_1}, readout_amp_2={self.readout_amp_2}, "
            f"RO_LO={self.RO_LO}, RO_LO_pwr={self.RO_LO_pwr}, ROq3={self.ROq3}, "
            f"ro_pulse_duration={self.ro_pulse_duration}, ro_dur={self.ro_dur})"
        )