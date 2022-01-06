import base.events
import user_sites


# TODO define single and double adsorption - ask if triple is pertinent?
class COAdsorption(base.events.EventBase):
    """
    Assume CO adsorbs on a top site
    """
    def __init__(self, params):
        base.events.EventBase.__init__(self, params, name='COAdsorption')

    def possible(self, system, site, other_site):
        # The site must be vacant and it must be a top site
        if system.sites[site].covered == user_sites.SPECIES_X and \
                system.sites[site].stype == user_sites.SITE_FCC111_TOP:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1.0 * self.params['pCO']
        return R

    def do_event(self, system, site, other_site):
        # cover one sites with species 1 (dissociative adsorption)
        system.sites[site].covered = user_sites.SPECIES_COX

    def get_involve_other(self):
        return False


class O2Adsorption(base.events.EventBase):
    """
    # TODO figure out which sites O2 desorbs on
    """
    def __init__(self, params):
        base.events.EventBase.__init__(self, params, name='O2Adsorption')

    def possible(self, system, site, other_site):
        if system.sites[site].covered == user_sites.SPECIES_X and \
                system.sites[other_site].covered == user_sites.SPECIES_X:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1.0 * self.params['pO2']
        return R

    def do_event(self, system, site, other_site):
        system.sites[site].covered = user_sites.SPECIES_OX
        system.sites[other_site].covered = user_sites.SPECIES_OX

    def get_involve_other(self):
        return False


class CODesorption(base.events.EventBase):
    def __init__(self, params):
        base.events.EventBase.__init__(self, params, name='CODesorption')

    def possible(self, system, site, other_site):
        if system.sites[site].covered == user_sites.SPECIES_COX:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1.0
        return R

    def do_event(self, system, site, other_site):
        system.sites[site].covered = user_sites.SPECIES_X

    def get_involve_other(self):
        return False


class O2Desorption(base.events.EventBase):
    def __init__(self, params):
        base.events.EventBase.__init__(self, params, name='O2Desorption')

    def possible(self, system, site, other_site):
        if system.sites[site].covered == user_sites.SPECIES_OX and \
                system.sites[other_site].covered == user_sites.SPECIES_OX:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1.0
        return R

    def do_event(self, system, site, other_site):
        system.sites[site].covered = user_sites.SPECIES_X
        system.sites[other_site].covered = user_sites.SPECIES_X

    def get_involve_other(self):
        return False


class CODiffusion(base.events.EventBase):
    def __init__(self, params):
        base.events.EventBase.__init__(self, params, name='CODiffusion')

    def possible(self, system, site, other_site):
        if system.sites[site].covered == user_sites.SPECIES_COX and \
                system.sites[other_site].covered == user_sites.SPECIES_X:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1.0
        return R
        # TODO compute diffusion rate

    def do_event(self, system, site, other_site):
        system.sites[other_site].covered = user_sites.SPECIES_COX
        system.sites[site].covered = user_sites.SPECIES_X

    def get_involve_other(self):
        return False


class ODiffusion(base.events.EventBase):
    def __init__(self, params):
        base.events.EventBase.__init__(self, params, name='ODiffusion')

    def possible(self, system, site, other_site):
        if system.sites[site].covered == user_sites.SPECIES_OX and \
                system.sites[other_site].covered == user_sites.SPECIES_X:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1.0
        return R
        # TODO get diffusion rate

    def do_event(self, system, site, other_site):
        system.sites[other_site].covered = user_sites.SPECIES_OX
        system.sites[site].covered = user_sites.SPECIES_X

    def get_involve_other(self):
        return False
