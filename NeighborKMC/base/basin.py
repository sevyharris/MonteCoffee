"""Methods for performing temporal acceleration of kMC simulations.

The methods here aid in performing the acceleration of the kMC simulations in MonteCoffee.
This is mainly based on the work of Dybeck et al. (https://doi.org/10.1021/acs.jctc.6b00859)

"""

import numpy as np


def rescaling(sim):
    """Rescales the times of occurrences for events.

    Rescales the times according to each quasi-equilibrated
    events *alpha*.

    Parameters
    ----------
    sim: NeighborKMC object
        main simulator object to perform rescaling of events for.

    """
    for ev in sim.equilEV:
        # Raise the barrier
        i_up = [i for i in range(len(sim.evs)) if sim.evs[i] == ev]
        for i in i_up:
            site = sim.siteslist[i]  # The site to do event.
            othersite = sim.other_sitelist[i]
            poss = sim.events[sim.evs[i]].possible(sim.system,
                                                     site, othersite)
            if poss:
                sim.rs[i] = sim.events[sim.evs[i]]. \
                    get_rate(sim.system, site, othersite)
                try:
                    u0 = -np.log(sim.us[i]) / sim.rs[i]
                    if sim.t < sim.tgen[i] + u0:
                        sim.frm_times[i] = sim.tgen[i] + u0
                    else:
                        sim.us[i] = uniform(0, 1)
                        sim.tgen[i] = sim.t
                        sim.frm_times[i] = sim.t - \
                                            np.log(sim.us[i]) / sim.rs[i]
                except:
                    sim.frm_times[i] = sim.tinfinity


def leave_superbasin(sim):
    """Leaves the superbasin.

    Resets all rate-scalings and statistics
    connected to the superbasin.

    Parameters
    ----------
    sim: NeighborKMC object
        main simulator object to perform rescaling of events for.

    """

    for e in sim.equilEV:
        sim.events[e].alpha = 1.
    rescaling(sim)

    sim.Suffex = []
    sim.r_S = np.zeros(len(sim.events))
    sim.dt_S = []
    sim.nem = np.zeros(len(sim.events), dtype=int)
    sim.isup = 0


def scaling_ks(sim, noneqevents, dtS):
    """Computes rate-constant based superbasin escape rate.

    Calculates superbasin escape rate according to the maximal rate-constant of
    events escaping the superbasin. This can be good for stability of the time-step.

    Parameters
    ----------
    sim: NeighborKMC object
        main simulator object to perform rescaling of events for.

    noneqevents: array-like
        The indices of events that are not in equilibrium, according to the loading
        order passed to *sim*.

    dtS: float
        The time spent in current superbasin (not used, but here for technical reasons as
        scaling_rs() needs this argument.)

    Returns
    --------
    float
        The estimated superbasin escape-rate.

    """
    if len(noneqevents)>0:
        r_S = max([sim.ksavg[neqev] for neqev in noneqevents])
    else:  # If all events are equilibrated, scale based on all rate consts.
        r_S = max(sim.ksavg)
    return r_S


def scaling_rs(sim, noneqevents, dtS):
    """Rate based superbasin escape time.

    Calculates superbasin escape time
    according to non-equilibrated event rates escaping
    the superbasin.

    c.f. The generalized temporal acceleration scheme
    (DOI: 10.1021/acs.jctc.6b00859)

    Parameters
    ----------
    sim: NeighborKMC object
        main simulator object to perform rescaling of events for.

    noneqevents: array-like
        The indices of events that are not in equilibrium, according to the loading
        order passed to *sim*.

    dtS: float
        The time spent in current superbasin (not used, but here for technical reasons as
        scaling_rs() needs this argument.)

    Returns
    --------
    float
        The estimated superbasin escape-rate.

    """
    if len(noneqevents) > 0:
        r_S = sum([sim.r_S[neqev] for neqev in noneqevents]) / dtS
    else:
        r_S = max(sim.r_S)

    return r_S


def superbasin(sim, evtype, dt):
    """Scales rates or leaves the current superbasin.

    Based on the current Monte Carlo step, the method determines
    if the superbasin is left or the fast quasi-equilibrated events
    should be slowed down.

    Keeps track and performs barrier adjustments,
    of the generalized temporal acceleration scheme
    (DOI: 10.1021/acs.jctc.6b00859)

    Parameters
    ----------
    sim: NeighborKMC object
        main simulator object to perform rescaling of events for.

    evtype: int
        The index of the event-type of the currently attempted Monte Carlo step.

    dt: float
        The time-step increment of a result of the currently attempted Monte Carlo step.


    Raises
    --------
    Warning
        If the time-step dt is negative.
    """
    # Update the rates in the current superbasin
    if dt < 0:
        raise Warning("Time-step is < 0. Are the events and neighborlists correct?. Exiting!!")
    farg = int(sim.frm_arg)
    sim.pm = (sim.pm + 1) % sim.ne
    sim.nem[evtype] += 1.
    sim.Nm[evtype][sim.pm] = 1.

    sim.r_S += [(sim.rs * dt)[sim.wheres[i][0]].sum() for i in range(len(sim.events))]
    sim.dt_S.append(dt)

    # See if event is quasi-equilibrated
    if evtype in sim.reverses:
        rev = abs(sim.Nm[evtype].sum() - \
                  sim.Nm[sim.reverses[evtype]].sum())

        Nexm = sim.Nm[evtype].sum() + \
               sim.Nm[sim.reverses[evtype]].sum()

        if evtype not in sim.equilEV:
            if Nexm >= sim.ne / 2. and rev < sim.delta * sim.ne:
                sim.equilEV.append(evtype)
                if evtype != sim.reverses[evtype]:
                    sim.equilEV.append(sim.reverses[evtype])

            else:
                leave_superbasin(sim)

        if evtype in sim.equilEV and sim.nem[evtype] + \
                sim.nem[sim.reverses[evtype]] >= sim.ne \
                and evtype in sim.equilEV \
                and evtype not in sim.Suffex:

            sim.Suffex.append(evtype)

            if evtype != sim.reverses[evtype]:
                sim.Suffex.append(sim.reverses[evtype])

    else:  # Not reversible
        leave_superbasin(sim)

    if sim.isup > sim.Ns:  # If observation period is over, scale events.
        dtS = sum(sim.dt_S)
        E = [i for i in range(len(sim.events)) if i not in sim.Suffex]
        r_S = sim.scaling_func(sim, E, dtS)

        for ev in [e for e in sim.equilEV if e in sim.Suffex]:
            rmev = sim.r_S[ev] / dtS
            rmrev = sim.r_S[sim.reverses[ev]] / dtS

            alpham = min(sim.Nf * r_S / (rmev + rmrev), 1)
            sim.events[ev].alpha *= alpham

        rescaling(sim)

        sim.isup = 0

    sim.isup += 1
