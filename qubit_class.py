class Qubit:
    def __init__(self, config):
        """
        Initialize a Qubit instance with specific properties from a dictionary.

        Args:
            config (dict): A dictionary containing key-value pairs for the Qubit properties.

        Raises:
            ValueError: If a key in the dictionary is not an allowed attribute.
        """
        # Define the allowed attributes for a Qubit
        self.allowed_attributes = {
            "qubit_id",
            "ro_freq",
            "ge_time",
            "ef_time",
            "ge_ssm",
            "ef_ssm",
            "ge_amp",
        }

        # Validate and set attributes from the dictionary
        for key, value in config.items():
            if key not in self.allowed_attributes:
                raise ValueError(
                    f"Invalid attribute '{key}' for Qubit. Allowed attributes are: {self.allowed_attributes}"
                )
            setattr(self, key, value)

        # Set default values (None) for missing attributes
        for attr in self.allowed_attributes:
            if not hasattr(self, attr):
                setattr(self, attr, None)

    def __repr__(self):
        # Dynamically create a string representation based on instance attributes
        attributes = ", ".join(
            f"{attr}={getattr(self, attr)}" for attr in self.allowed_attributes
        )
        return f"Qubit({attributes})"


class Readout:
    def __init__(self, config):
        """
        Initialize a Readout instance with specific properties from a dictionary.

        Args:
            config (dict): A dictionary containing key-value pairs for the Readout properties.

        Raises:
            ValueError: If a key in the dictionary is not an allowed attribute.
        """
        # Define the allowed attributes for a Readout
        self.allowed_attributes = {
            "readout_amp_1",
            "readout_amp_2",
            "RO_LO",
            "RO_LO_pwr",
            "ROq3",
            "ro_dur",
        }
        # Validate and set attributes from the dictionary
        for key, value in config.items():
            if key not in self.allowed_attributes:
                raise ValueError(
                    f"Invalid attribute '{key}' for Readout. Allowed attributes are: {self.allowed_attributes}"
                )
            setattr(self, key, value)

        # Set default values (None) for missing attributes
        for attr in self.allowed_attributes:
            if not hasattr(self, attr):
                setattr(self, attr, None)

    def __repr__(self):
        # Dynamically create a string representation based on instance attributes
        attributes = ", ".join(
            f"{attr}={getattr(self, attr)}" for attr in self.allowed_attributes
        )
        return f"Readout({attributes})"