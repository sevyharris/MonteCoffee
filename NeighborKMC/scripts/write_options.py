import ConfigParser

config = ConfigParser.RawConfigParser()
config.add_section('Parameters')
config.set('Parameters', 'Nspecies', '2')
config.set('Parameters', 'NNinteractions', '1')
config.set('Parameters', 'SaveSteps', '100000') # Save a pickle file every 100000 steps.
config.set('Parameters', 'LogSteps', '1000')
config.set('Parameters', 'PicklePrefix', 'transitions_step')

 

config.add_section('Options')
config.set('Options', 'SaveCovs', True)
config.set('Options', 'Verbose', True)

# Writing our configuration file to 'example.cfg'
with open('kMC_options.cfg', 'wb') as configfile:
    config.write(configfile)
