"""Implement OperatingConditions class.

This module contains the OperatingConditions class used to
store information regarding the operating conditions
in freezing processes.
"""
import numpy as np
import warnings

from typing import Optional, Iterable, Union


class OperatingConditions:
    """A class to handle a single Stochastic Nucleation of Water simulation.

    More information regarding the equations and their derivation can be found in
    "Stochastic shelf-scale modeling framework for the freezing stage in
    freeze-drying processes", Deck, Ochsenbein, and Mazzotti (2022),
    Int J Pharm, 613, 121276, https://doi.org/10.1016/j.ijpharm.2021.121276.

    Parameters:
        cnt (float): Controlled nucleation time.
        controlledNucleation (bool): Controlled nucleation on/off.
        cooling (dict): A dictionary describing the cooling profile.
        holding (dict): A dictionary describing the holding step.
        t_tot (float): The total process time.
    """

    def __init__(
        self,
        t_tot: float = 2e4,
        cooling: dict = {"rate": 0.5 / 60, "start": 20, "end": -50},
        holding: Optional[Union[Iterable[dict], dict]] = None,
        cnTemp: Union[float, int] = None,
    ):
        """Construct an OperatingConditions object.

        Args:
            t_tot (float, optional): The total process time. Defaults to 2e4.
            cooling (dict, optional): A dictionary describing the cooling profile.
                Defaults to {"rate": 0.5 / 60, "start": 20, "end": -50}.
            holding (Optional[Union[Iterable[dict], dict]], optional):
                A dictionary or list of dictionaries describing
                the holding step(s). Defaults to None (no holding).
            cnTemp (Union[float, int], optional): At what temperature controlled
                nucleation is triggered. Nucleation triggers when
                temp<=cnTemp. If that occurs during a holding
                phase nucleation will trigger _at the end_ of the phase.
        """
        self.t_tot = t_tot
        if not all([key in cooling.keys() for key in ["rate", "start", "end"]]):
            raise ValueError("Cooling dictionary does not contain all required keys.")

        # deal with constant temp; basically a convenience wrapper
        if (cooling["rate"] == 0) and (cooling["start"] == cooling["end"]):
            cooling["rate"] = 1e-16  # rate equal zero would yields divide by zero error

            if holding is None:
                holding = dict(duration=t_tot, temp=cooling["start"])
            elif holding is not None:
                raise ValueError("Cannot have cooling rate 0 _and_ holding steps.")

        elif (cooling["rate"] == 0) and (cooling["start"] != cooling["end"]):
            raise ValueError(
                "If cooling rate is zero, start must match end temperature."
            )

        self.cooling = cooling

        self.holding = holding

        if self._t_tot_implied > t_tot:
            warnings.warn(
                (
                    f"The implied process time (inferred from cooling rate + holding step(s)) "
                    + f"is larger than the total time t_tot ({self._t_tot_implied}>{t_tot}). "
                    + f"The simulation will stop at t_tot={t_tot}."
                ),
                UserWarning,
            )

        self.cnTemp = cnTemp

    @property
    def holding(self) -> Iterable[dict]:
        """Get holding property."""
        return self._holding

    @property
    def _t_tot_implied(self) -> float:
        coolingTime = (self.cooling["start"] - self.cooling["end"]) / self.cooling[
            "rate"
        ]
        if (self.holding is not None) and (isinstance(self.holding, (tuple, list))):
            holdingTime = sum([h["duration"] for h in self.holding])
        elif (self.holding is not None) and (isinstance(self.holding, (dict))):
            holdingTime = self.holding["duration"]
        else:
            holdingTime = 0
        return coolingTime + holdingTime

    @holding.setter
    def holding(self, value: Union[dict, list, tuple]):
        """Set holding property."""
        if isinstance(value, dict):
            value = [value]
        elif isinstance(value, (list, tuple)):
            # bring holding steps into right order (descending)
            value = sorted(value, key=lambda hdict: hdict["temp"], reverse=True)
        elif value is None:
            pass
        else:
            raise TypeError("holding must be a dict or Iterable of dict.")

        if value is not None:
            for val in value:
                if not all([key in val.keys() for key in ["duration", "temp"]]):
                    raise ValueError(
                        "Holding dictionary does not contain all required keys."
                    )

        self._holding = value

    @property
    def cnt(self) -> float:
        """Return the time when controlled nucleation should trigger.

        Raises:
            NotImplementedError: If holding is not defined
                don't know how to calculate cnt.

        Returns:
            float: The time of controlled nucleation
                (inf if no controlled nucleation applied).
        """
        if self.cnTemp is not None:
            T_vec = self.tempProfile(1)
            t_vec = np.arange(0, len(T_vec))

            I_endHold = np.argmax(T_vec[::-1] >= self.cnTemp)
            t_endHold = t_vec[::-1][I_endHold]

        else:
            t_endHold = np.inf

        return t_endHold

    def tempProfile(self, dt: float) -> np.ndarray:
        """Return temperature profile.

        Compute temperature profile with or without
        holding step.
        Args:
            dt (float): The time step size.

        Returns:
            np.ndarray: The temperature profile.
        """
        # total number of steps
        n = int(np.ceil(self.t_tot / dt)) + 1

        T_start = self.cooling["start"]
        cr = self.cooling["rate"]

        if self.holding is not None:
            hdicts = self.holding + [
                {"temp": self.cooling["end"], "duration": self.t_tot}
            ]
        else:
            hdicts = [{"temp": self.cooling["end"], "duration": self.t_tot}]

        T_vec = np.array([])

        for hdict in hdicts:

            T_hold = hdict["temp"]
            t_hold = (T_start - T_hold) / cr
            duration_hold = hdict["duration"]

            # ramp from start to hold temp
            T_vec_toHold = self._simpleCool(
                Tstart=T_start,
                Tend=T_hold,
                coolingRate=cr,
                dt=dt,
            )

            # append holding period
            T_vec_holding = [T_hold] * int(np.ceil((duration_hold - t_hold % dt) / dt))

            T_start = T_hold

            T_vec = np.concatenate([T_vec, T_vec_toHold, T_vec_holding])

        T_vec = T_vec[:n]

        return T_vec

    def _simpleCool(
        self,
        Tstart: float,
        Tend: float,
        coolingRate: float,
        dt: float,
        t_tot: Optional[float] = None,
    ) -> np.ndarray:
        """Return cooling profile for linear step.

        Args:
            Tstart (float): Start temperature.
            Tend (float): End temperature.
            coolingRate (float): Cooling rate.
            dt (float): Time step.
            t_tot (Optional[float], optional): Total process time.
                Defaults to None.

        Returns:
            np.ndarray: The temperature profile.
        """
        t_end = (Tstart - Tend) / coolingRate
        t_vec = np.arange(0, t_end, dt)

        T_profile = Tstart - t_vec * coolingRate

        if t_tot is not None:
            assert (
                t_tot >= t_end
            ), "Final temp can't be reached with cooling rate, holding/total time."

            add_n = int(np.ceil((t_tot - t_vec[-1]) / dt))

            T_profile = np.append(T_profile, [Tend] * add_n)

        return T_profile

    def __repr__(self) -> str:
        """Return string representation of the OperatingConditions class.

        Returns:
            str: The OperatingConditions class string representation
                giving some basic info.
        """
        if self.holding is not None:
            holdPluralBool = len(self.holding) > 1
            holdPlural = f"Hold{'s' if (holdPluralBool) else ''}: "
            holdStr = holdPlural + " AND ".join(
                [f"{hdict['duration']} @ {hdict['temp']}" for hdict in self.holding]
            )
        else:
            holdStr = "No Holds"

        return (
            f"OperatingConditions([t_tot: {self.t_tot}, "
            + f"Cooling: {self.cooling['start']} to {self.cooling['end']} "
            + f"with rate {self.cooling['rate']:4.2f}, "
            + holdStr
            + ", "
            + f"Controlled Nucleation: {'ON' if self.cnTemp else 'OFF'}"
        )
