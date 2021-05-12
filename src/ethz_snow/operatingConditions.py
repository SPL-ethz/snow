'''
Created Date: Wednesday May 5th 2021
Author: David Ochsenbein (DRO) - dochsenb@its.jnj.com
-----
Copyright (c) 2021 David Ochsenbein, Johnson & Johnson
'''
import numpy as np

from typing import Optional


class OperatingConditions:
    def __init__(
            self,
            t_tot: float = 2e4,
            cooling: dict = {'rate': 0.5/60, 'start': 20, 'end': -50},
            holding: Optional[dict] = {'duration': 10, 'temp': -12},
            controlledNuclation: bool = False):

        self.t_tot = t_tot
        if not all([key in cooling.keys() for key in ['rate', 'start', 'end']]):
            raise ValueError("Cooling dictionary does not contain all required keys.")
        self.cooling = cooling

        if not all([key in holding.keys() for key in ['duration', 'temp']]):
            raise ValueError("Holding dictionary does not contain all required keys.")
        self.holding = holding

        self.controlledNuclation = controlledNuclation

    def tempProfile(self, dt: float) -> np.ndarray:

        # total number of steps
        n = int(np.ceil(self.t_tot/dt))+1

        if self.holding is not None:
            T_hold = self.holding['temp']
            duration_hold = self.holding['duration']

            # time and number of steps to hold temperature
            T_vec_toHold = self._simpleCool(Tstart=self.cooling['start'], Tend=T_hold,
                                            coolingRate=self.cooling['rate'], dt=dt)

            # append holding period
            T_vec_holding = [T_hold]*int(np.ceil(duration_hold / dt))

            T_vec_toEnd = self._simpleCool(Tstart=T_hold, Tend=self.cooling['end'],
                                           coolingRate=self.cooling['rate'],
                                           dt=dt, t_tot=self.t_tot)

            T_vec = np.concatenate([T_vec_toHold, T_vec_holding, T_vec_toEnd])

            T_vec = T_vec[:n]
        else:

            T_vec = self._simpleCool(Tstart=self.cooling['start'], Tend=self.cooling['end'],
                                     coolingRate=self.cooling['rate'], dt=dt, t_tot=self.t_tot)

        return T_vec

    def _simpleCool(
            self,
            Tstart: float,
            Tend: float,
            coolingRate: float,
            dt: float,
            t_tot: Optional[float] = None) -> np.ndarray:

        t_end = (Tstart - Tend) / coolingRate
        n_end = int(np.ceil(t_end / dt))

        t_vec = np.linspace(0, t_end, n_end)

        T_profile = Tstart - t_vec * coolingRate

        if t_tot is not None:
            assert t_tot >= t_end, \
                "Final Temperature cannot be reached with given cooling rate, holding/total time."

            add_n = int(np.ceil((t_tot - t_end) / dt))

            T_profile = np.append(T_profile, [Tend]*add_n)

        return T_profile

    def __repr__(self) -> str:
        """ The string representation of the OperatingConditions class.

        Returns:
            str: The OperatingConditions class string representation giving some basic info.
        """

        return (f"OperatingConditions([t_tot: {self.t_tot}, "
                + f"Cooling: {self.cooling['start']} to {self.cooling['end']} "
                + f"with {self.cooling['rate']}, "
                + f"Hold: {self.holding['duration']} @ {self.holding['temp']}, "
                + f"Controlled Nucleation: {'ON' if self.controlledNuclation else 'OFF'}")
