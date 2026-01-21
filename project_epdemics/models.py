from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class TwoStrainParams:
    """
    Parameter container for two-strain epidemic models.

    This class defines the biological and epidemiological parameters
    shared across deterministic (ODE) and stochastic formulations.

    Parameters
    ----------
    beta1, beta2 : float
        Transmission rates for strain 1 and strain 2.

    gamma1, gamma2 : float
        Recovery rates for strain 1 and strain 2.

    omega1, omega2 : float, optional
        Immunity waning rates (SIRS). Set to 0 for SIR models.

    sigma12 : float, optional
        Cross-immunity parameter describing susceptibility to strain 2
        after recovery from strain 1.
        sigma12 = 0   → complete cross-protection
        sigma12 = 1   → no cross-protection

    sigma21 : float, optional
        Cross-immunity parameter describing susceptibility to strain 1
        after recovery from strain 2.

    N : float
        Total population size (assumed constant).
    """
    beta1: float
    beta2: float
    gamma1: float
    gamma2: float
    omega1: float = 0.0
    omega2: float = 0.0
    sigma12: float = 1.0
    sigma21: float = 1.0
    N: float = 1.0


def rhs_two_strain_sirs(t: float, y: np.ndarray, p: TwoStrainParams) -> np.ndarray:
    """
    Right-hand side of a two-strain SIRS model.

    Compartments
    ------------
    S  : Susceptible individuals
    I1 : Individuals infected with strain 1
    I2 : Individuals infected with strain 2
    R1 : Recovered from strain 1
    R2 : Recovered from strain 2

    Modeling assumptions
    --------------------
    - Well-mixed population
    - No co-infection
    - Partial cross-immunity implemented as reduced susceptibility
    - Deterministic mean-field approximation of a stochastic process
    """

    S, I1, I2, R1, R2 = y
    N = p.N

    # Forces of infection (per-susceptible infection pressure)
    lam1 = p.beta1 * I1 / N
    lam2 = p.beta2 * I2 / N

    # Primary infections from the susceptible pool
    inf1_S = lam1 * S
    inf2_S = lam2 * S

    # Cross-infections from recovered states (partial cross-immunity)
    inf2_R1 = p.sigma12 * lam2 * R1
    inf1_R2 = p.sigma21 * lam1 * R2

    # Recovery processes
    rec1 = p.gamma1 * I1
    rec2 = p.gamma2 * I2

    # Immunity waning (SIRS only)
    wan1 = p.omega1 * R1
    wan2 = p.omega2 * R2

    # Differential equations
    dS  = -inf1_S - inf2_S + wan1 + wan2
    dI1 = +inf1_S + inf1_R2 - rec1
    dI2 = +inf2_S + inf2_R1 - rec2
    dR1 = +rec1 - wan1 - inf2_R1
    dR2 = +rec2 - wan2 - inf1_R2

    return np.array([dS, dI1, dI2, dR1, dR2], dtype=float)


def rhs_two_strain_sis(t: float, y: np.ndarray, p: TwoStrainParams) -> np.ndarray:
    """
    Right-hand side of a two-strain SIS model.

    Compartments
    ------------
    S  : Susceptible individuals
    I1 : Individuals infected with strain 1
    I2 : Individuals infected with strain 2

    Notes
    -----
    - Recovery returns individuals directly to the susceptible class
    - No explicit immune memory or cross-immunity
    - Useful as a simplified baseline model
    """

    S, I1, I2 = y
    N = p.N

    lam1 = p.beta1 * I1 / N
    lam2 = p.beta2 * I2 / N

    inf1 = lam1 * S
    inf2 = lam2 * S

    rec1 = p.gamma1 * I1
    rec2 = p.gamma2 * I2

    dS  = -inf1 - inf2 + rec1 + rec2
    dI1 = +inf1 - rec1
    dI2 = +inf2 - rec2

    return np.array([dS, dI1, dI2], dtype=float)


def check_nonnegative(y: np.ndarray, atol: float = 1e-9) -> bool:
    """
    Numerical sanity check ensuring all compartments remain non-negative.
    Useful for detecting solver or modeling errors.
    """
    return bool(np.all(y >= -atol))


def check_mass_conservation(y: np.ndarray, N: float, atol: float = 1e-6) -> bool:
    """
    Check conservation of total population mass.
    This should hold exactly for closed epidemic models.
    """
    return bool(abs(np.sum(y) - N) <= atol)