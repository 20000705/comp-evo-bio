from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from numpy.random import Generator, default_rng
from scipy.integrate import solve_ivp

from models import TwoStrainParams, check_nonnegative, check_mass_conservation


# ============================================================
# Deterministic ODE solver wrapper
# ============================================================

@dataclass
class ODEResult:
    """
    Container for deterministic ODE results.

    Attributes
    ----------
    t : np.ndarray
        Time points returned by the solver.

    y : np.ndarray
        State trajectories with shape (dim, len(t)).
        Each row corresponds to a compartment.
    """
    t: np.ndarray
    y: np.ndarray


def run_ode(
    rhs: Callable[[float, np.ndarray, TwoStrainParams], np.ndarray],
    y0: np.ndarray,
    t_span: Tuple[float, float],
    p: TwoStrainParams,
    t_eval: Optional[np.ndarray] = None,
    rtol: float = 1e-6,
    atol: float = 1e-9,
) -> ODEResult:
    """
    Solve a deterministic two-strain epidemic model using SciPy's ODE solver.

    This function represents the *mean-field limit* of the epidemic process.
    It is used to study average epidemic dynamics and equilibria.

    Parameters
    ----------
    rhs : callable
        Right-hand side function of the ODE system.

    y0 : np.ndarray
        Initial condition (compartment densities or counts).

    t_span : (float, float)
        Time interval (t0, tf).

    p : TwoStrainParams
        Model parameters.

    t_eval : np.ndarray, optional
        Specific times at which to store the solution.

    Returns
    -------
    ODEResult
        Time points and solution trajectories.
    """
    sol = solve_ivp(
        fun=lambda t, y: rhs(t, y, p),
        t_span=t_span,
        y0=y0.astype(float),
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
    )

    if not sol.success:
        raise RuntimeError(f"ODE solve failed: {sol.message}")

    return ODEResult(t=sol.t, y=sol.y)


# ============================================================
# Stochastic Gillespie simulation
# ============================================================

@dataclass
class GillespieEvent:
    """
    Record of a single stochastic event.

    Attributes
    ----------
    t : float
        Time at which the event occurred.

    event : str
        Event type (e.g., 'S->I1', 'I1->R1').

    state : np.ndarray
        System state immediately after the event.
    """
    t: float
    event: str
    state: np.ndarray


@dataclass
class GillespieResult:
    """
    Container for Gillespie simulation output.

    Attributes
    ----------
    times : np.ndarray
        Recorded event times (not uniformly spaced).

    states : np.ndarray
        Recorded states with shape (n_records, dim).

    events : list[GillespieEvent]
        Full event log.

    tree_parents, tree_times, tree_labels :
        Placeholders for future transmission / evolutionary tree construction.
    """
    times: np.ndarray
    states: np.ndarray
    events: List[GillespieEvent]

    tree_parents: Optional[List[int]] = None
    tree_times: Optional[List[float]] = None
    tree_labels: Optional[List[str]] = None


def run_gillespie_two_strain_sirs(
    y0: np.ndarray,
    t_max: float,
    p: TwoStrainParams,
    rng: Optional[Generator] = None,
    max_steps: int = 2_000_000,
    record_every: int = 1,
    sanity_checks: bool = True,
) -> GillespieResult:
    """
    Gillespie stochastic simulation algorithm (SSA) for a two-strain SIRS model.

    This model explicitly tracks *integer counts* of individuals and captures
    stochastic effects such as extinction and lineage loss, which are not
    represented in deterministic ODE models.

    State vector
    ------------
    y = [S, I1, I2, R1, R2]

    Event types
    -----------
    Primary infection:
      S -> I1
      S -> I2

    Recovery:
      I1 -> R1
      I2 -> R2

    Waning immunity:
      R1 -> S
      R2 -> S

    Cross-infection (partial cross-immunity):
      R1 -> I2
      R2 -> I1

    Notes
    -----
    - The stochastic model converges to the ODE solution in expectation as N grows.
    - Stochasticity is essential for studying lineage persistence and tree structure.
    """

    rng = default_rng() if rng is None else rng

    # Initialize integer state
    y = y0.astype(int).copy()
    if y.shape[0] != 5:
        raise ValueError("Expected y0 = [S, I1, I2, R1, R2].")

    # Consistency checks
    if sanity_checks:
        if int(np.sum(y)) != int(p.N):
            raise ValueError(
                f"Population mismatch: sum(y0)={np.sum(y)} but p.N={p.N}"
            )
        if not check_nonnegative(y):
            raise ValueError("Initial state contains negative values.")

    # Initialize time and storage
    t = 0.0
    times = [t]
    states = [y.copy()]
    events: List[GillespieEvent] = []

    # --------------------------------------------------------
    # Event rate (hazard) computation
    # --------------------------------------------------------
    def rates(y: np.ndarray) -> Tuple[np.ndarray, List[str]]:
        """
        Compute event rates given the current state.

        Returns
        -------
        r : np.ndarray
            Event rate vector.

        names : list[str]
            Names corresponding to each rate.
        """
        S, I1, I2, R1, R2 = y
        N = p.N

        # Force of infection
        lam1 = p.beta1 * I1 / N
        lam2 = p.beta2 * I2 / N

        r = np.array([
            lam1 * S,              # S -> I1
            lam2 * S,              # S -> I2
            p.gamma1 * I1,         # I1 -> R1
            p.gamma2 * I2,         # I2 -> R2
            p.omega1 * R1,         # R1 -> S
            p.omega2 * R2,         # R2 -> S
            p.sigma12 * lam2 * R1, # R1 -> I2
            p.sigma21 * lam1 * R2, # R2 -> I1
        ], dtype=float)

        names = [
            "S->I1", "S->I2",
            "I1->R1", "I2->R2",
            "R1->S", "R2->S",
            "R1->I2", "R2->I1",
        ]

        return r, names

    # --------------------------------------------------------
    # State update for each event
    # --------------------------------------------------------
    def apply_event(y: np.ndarray, ev: str) -> None:
        """Apply one event to the state vector in-place."""
        if ev == "S->I1":
            y[0] -= 1; y[1] += 1
        elif ev == "S->I2":
            y[0] -= 1; y[2] += 1
        elif ev == "I1->R1":
            y[1] -= 1; y[3] += 1
        elif ev == "I2->R2":
            y[2] -= 1; y[4] += 1
        elif ev == "R1->S":
            y[3] -= 1; y[0] += 1
        elif ev == "R2->S":
            y[4] -= 1; y[0] += 1
        elif ev == "R1->I2":
            y[3] -= 1; y[2] += 1
        elif ev == "R2->I1":
            y[4] -= 1; y[1] += 1
        else:
            raise ValueError(f"Unknown event: {ev}")

    # --------------------------------------------------------
    # Main Gillespie loop
    # --------------------------------------------------------
    for step in range(max_steps):
        if t >= t_max:
            break

        r, names = rates(y)
        r_tot = float(np.sum(r))
        if r_tot <= 0.0:
            break  # epidemic extinct or frozen

        # Sample waiting time
        dt = rng.exponential(1.0 / r_tot)
        t_next = t + dt
        if t_next > t_max:
            t = t_max
            break

        # Sample event
        cdf = np.cumsum(r)
        u = rng.random() * r_tot
        idx = int(np.searchsorted(cdf, u, side="right"))
        if idx >= len(names):  # numerical edge case
            idx = len(names) - 1

        ev = names[idx]

        # Apply event
        apply_event(y, ev)
        t = t_next

        # Safety checks
        if sanity_checks:
            if not check_nonnegative(y):
                raise RuntimeError(f"Negative state after {ev} at t={t}")
            if not check_mass_conservation(y, p.N):
                raise RuntimeError(f"Mass not conserved at t={t}")

        # Record trajectory (event-based sampling)
        if (step % record_every) == 0:
            times.append(t)
            states.append(y.copy())

        events.append(GillespieEvent(t=t, event=ev, state=y.copy()))

    return GillespieResult(
        times=np.array(times, dtype=float),
        states=np.vstack(states).astype(int),
        events=events,
    )