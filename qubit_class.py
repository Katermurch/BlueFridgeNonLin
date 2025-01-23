class Qubit:
    def __init__(
        self,
        qubit_id,
        ro_freq=None,
        ge_time=None,
        ef_time=None,
        ge_ssm=None,
        ef_ssm=None,
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
        self.qubit_id = qubit_id
        self.ro_freq = ro_freq
        self.ge_time = ge_time
        self.ef_time = ef_time
        self.ge_ssm = ge_ssm
        self.ef_ssm = ef_ssm

    def __repr__(self):
        return (
            f"Qubit({self.qubit_id}, ro_freq={self.ro_freq}, ge_time={self.ge_time}, "
            f"ef_time={self.ef_time}, ge_ssm={self.ge_ssm}, ef_ssm={self.ef_ssm})"
        )
