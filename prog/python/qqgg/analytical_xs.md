
# Table of Contents

1.  [Init](#org6d91f18)
    1.  [Required Modules](#org84213bf)
    2.  [Utilities](#orgffa1831)
2.  [Implementation](#orge6ee38e)
3.  [Calculations](#org3896b17)
    1.  [XS qq -> gamma gamma](#orgaad68f9)


<a id="org6d91f18"></a>

# Init


<a id="org84213bf"></a>

## Required Modules

    import numpy as np
    import matplotlib.pyplot as plt

    [....]


<a id="orgffa1831"></a>

## Utilities

    %run ../utility.py

[&#x2026;.]


<a id="orge6ee38e"></a>

# Implementation

    """
    Implementation of the analytical cross section for q q_bar ->
    gamma gamma
    
    Author: Valentin Boettcher <hiro@protagon.space>
    """
    
    import numpy as np
    from scipy.constants import alpha
    
    # NOTE: a more elegant solution would be a decorator
    def energy_factor(charge, esp):
        """
        Calculates the factor common to all other values in this module
    
        Arguments:
        esp -- center of momentum energy in GeV
        charge -- charge of the particle in units of the elementary charge
        """
    
        return charge**4*(alpha/esp)**2/6
    
    
    def diff_xs(theta, charge, esp):
        """
        Calculates the differential cross section as a function of the
        azimuth angle theta in units of 1/GeV^2.
    
        Arguments:
        theta -- azimuth angle
        esp -- center of momentum energy in GeV
        charge -- charge of the particle in units of the elementary charge
        """
    
        f = energy_factor(charge, esp)
        return f*(1 + 2/np.sin(theta)**2)
    
    def diff_xs_eta(eta, charge, esp):
        """
        Calculates the differential cross section as a function of the
        pseudo rapidity of the photons in units of 1/GeV^2.
    
        Arguments:
        eta -- pseudo rapidity
        esp -- center of momentum energy in GeV
        charge -- charge of the particle in units of the elementary charge
        """
    
        f = energy_factor(charge, esp)
        return f*(1 + 2*np.cosh(eta)**2)
    
    def total_xs_eta(eta, charge, esp):
        """
        Calculates the total cross section as a function of the pseudo
        rapidity of the photons in units of 1/GeV^2.  If the rapditiy is
        specified as a tuple, it is interpreted as an interval.  Otherwise
        the interval [-eta, eta] will be used.
    
        Arguments:
        eta -- pseudo rapidity (tuple or number)
        esp -- center of momentum energy in GeV
        charge -- charge of the particle in units of the elementar charge
        """
    
        f = energy_factor(charge, esp)
        if not isinstance(eta, tuple):
            eta = (-eta, eta)
    
        if len(eta) != 2:
            raise ValueError('Invalid eta cut.')
    
        def F(x):
            return np.tanh(x) - 2*x
    
        return 2*np.pi*f*(F(eta[0]) - F(eta[1]))


<a id="org3896b17"></a>

# Calculations


<a id="orgaad68f9"></a>

## XS qq -> gamma gamma

First, set up the input parameters.

    eta = 2.5
    charge = 1/3
    esp = 200  # GeV

[&#x2026;.]

And now calculate the cross section in picobarn.

    xs_gev = total_xs_eta(eta, charge, esp)
    xs_pb = gev_to_pb(xs_gev)
    xs_pb

[&#x2026;.]

Compared to sherpa, its a little off.

    sherpa = 0.0538009
    xs_pb/sherpa

[&#x2026;.]

