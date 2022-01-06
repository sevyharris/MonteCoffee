"""Contains methods for performing temporal acceleration of kMC simulations.

The methods here aid in performing the acceleration of the kMC simulations in MonteCoffee.
This is mainly based on the work of Dybeck et al. (https://doi.org/10.1021/acs.jctc.6b00859)

"""

import numpy as np
from statistics import mean


def rescaling(sim):
    """Rescales the times of occurrences for events.

    Rescales the times according to each quasi-equilibrated and sufficiently executed
    events *alpha*.

    Parameters
    ----------
    sim: NeighborKMC
        main simulator object to perform rescaling of events for.

    """
    if len(sim.Suffex) > 0:
        for ev in sim.Suffex:
            # Raise the barrier
            i_up = [iu for iu in range(len(sim.evs)) if sim.evs[iu] == ev and sim.possible_evs[iu] and not iu in sim.executed_poslist]
            for i in i_up:
                site = sim.siteslist[i]  # The site to do event.
                othersite = sim.other_sitelist[i]
                old_rs = sim.rs[i]
                sim.rs[i] = sim.events[sim.evs[i]].get_rate(sim.system, site, othersite)
                try:
                    u0 = -np.log(sim.us[i]) / sim.rs[i]
                    if sim.t < sim.tgen[i] + u0:
                        sim.frm_times[i] = sim.tgen[i] + u0
                    else:
                        sim.us[i] = random.uniform(0, 1)
                        sim.tgen[i] = sim.t
                        sim.frm_times[i] = sim.t - np.log(sim.us[i]) / sim.rs[i]
                except:
                    sim.frm_times[i] = sim.tinfinity

        if sim.events[sim.Suffex[0]].alpha == 1:
            sim.Suffex = []


def leave_superbasin(sim):
    """Leaves the superbasin.

    Resets all rate-scalings and statistics
    connected to the superbasin. The sufficiently executed event list is
    reset with rescaling.

    Parameters
    ----------
    sim: NeighborKMC
        main simulator object to perform rescaling of events for.

    """

    for e in sim.equilEV:
        sim.events[e].alpha = 1.

    sim.r_S = np.zeros(len(sim.events))
    sim.k_S = np.zeros(len(sim.events))
    sim.dt_S = []
    sim.nem = np.zeros(len(sim.events), dtype=int)
    sim.isup = 0
    sim.equilEV = []


def scale_rate(sim):
    """Rate based superbasin escape time.

    Calculates superbasin escape time
    according to non-equilibrated event rates escaping
    the superbasin.

    c.f. the generalized temporal acceleration scheme
    of Dybeck et al.

    Parameters
    ----------
    sim: NeighborKMC
        main simulator object to perform rescaling of events for.

    noneqevents: list(int)
        The indices of events that are not in equilibrium, according to the loading
        order passed to *sim*.

    Returns
    --------
    float
        The estimated superbasin escape-rate.

    """
    do_scaling = False
    noneqevents = [i for i in range(len(sim.events)) if i not in sim.Suffex]
    if len(noneqevents) > 0:
        r_S = sum([sim.r_S[neqev] for neqev in noneqevents])
    else:
        r_S = max(sim.r_S)

    for ev in [e for e in sim.equilEV if e in sim.Suffex]:
        rmev = sim.r_S[ev]
        rmrev = sim.r_S[sim.reverses[ev]]
        alpham = min(2.0 * r_S / (rmev + rmrev) * sim.Nf, 1)
        print(ev, alpham, 2.0 * r_S / (rmev + rmrev))
        sim.events[ev].alpha = alpham
        do_scaling = True
    return do_scaling


def scale_constant(sim):
    """Rate based on a constant frequency factor for deeceleration
    as used for example by Hoffmann and Bligaard

    """
    do_scaling = False
    for ev in [e for e in sim.equilEV if e in sim.Suffex]:
        alpham = min(1.0 / float(sim.Nf), 1.)
        sim.events[ev].alpha = alpham
        do_scaling = True
    return do_scaling


def scale_rate_constant(sim):
    """Rates based on the mean of the current observation period """

    do_scaling = False
    noneqevents = [i for i in range(len(sim.events)) if i not in sim.Suffex and not sim.nem[i] == 0]
    if len(noneqevents) > 0:
        k_S = sum([(float(sim.k_S[neqev]) / float(sum(sim.nem))) for neqev in noneqevents])
    else:
        k_S = max(sim.k_S) / float(sum(sim.nem))

    for ev in [e for e in sim.equilEV if e in sim.Suffex]:
        kmev = (sim.k_S[ev] / sum(sim.nem))
        kmrev = (sim.k_S[sim.reverses[ev]] / sum(sim.nem))
        alpham = min(2.0 * k_S / (kmev + kmrev) * sim.Nf, 1)
        # print ( 2.*k_S /(kmev+ kmrev),  alpham,ev)
        sim.events[ev].alpha = alpham
        do_scaling = True
    return do_scaling


def superbasin(sim, evtype, dt):
    """Scales rates or leaves the current superbasin.

    Based on the current Monte Carlo step, the method determines
    if the superbasin is left or the fast quasi-equilibrated events
    should be slowed down.

    Keeps track and performs barrier adjustments,
    of the generalized temporal acceleration scheme of Dybeck et al.

    Parameters
    ----------
    sim: NeighborKMC
        main simulator object to perform rescaling of events for.

    evtype: int
        The index of the event-type of the currently attempted Monte Carlo step.

    dt: float
        The time-step of the currently attempted Monte Carlo step.

    do_scaling: bool
        Bool if the rate constants are to be rescaled or not.


    Raises
    --------
    Warning
        If the time-step is negative.
    """
    # Update the rates in the current superbasin
    if dt < 0:
        raise Warning("Time-step is < 0. Are the events and neighborlists correct?. Exiting!!")

    sim.nem[evtype] += 1.

    sim.r_S += [(sim.rs * dt)[sim.wheres[i][0]].sum() for i in range(len(sim.events))]
    sim.k_S += [(sim.rs)[sim.wheres[i][0]].mean() for i in range(len(sim.events))]

    do_scaling = False

    # See if event is quasi-equilibrated
    if evtype in sim.reverses:

        rev = abs(sim.nem[evtype] - sim.nem[sim.reverses[evtype]])/(sim.nem[evtype] + sim.nem[sim.reverses[evtype]])

        if evtype not in sim.equilEV:
            if rev < sim.delta:
                sim.equilEV.append(evtype)
                if evtype != sim.reverses[evtype]:
                    sim.equilEV.append(sim.reverses[evtype])

        if evtype in sim.equilEV and sim.nem[evtype] >= sim.ne and sim.nem[sim.reverses[evtype]] >= sim.ne \
                and evtype not in sim.Suffex:

            sim.Suffex.append(evtype)

            if evtype != sim.reverses[evtype]:
                sim.Suffex.append(sim.reverses[evtype])

    else:  # Not reversible
        leave_superbasin(sim)
        do_scaling = True

    if sim.isup > sim.Ns:  # If observation period is over, scale events.
        possibles = globals().copy()
        possibles.update(locals())
        do_scaling = possibles.get(sim.use_scaling_algorithm)(sim)
        sim.isup = 0

    sim.isup += 1
    return do_scaling
