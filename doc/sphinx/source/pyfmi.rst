

================================
Assimulo for PyFMI users
================================


Assimulo has been incorporated into PyFMI as the default simulation package. This has been made possible by extending Assimulo's problem classes where the FMU model is modified and adapted to Assimulo.

When using Assimulo together with PyFMI and the new high-level simulation methods, all the parameters passed in the *options*-dict represents parameters for the specific solver used in Assimulo. ::
 
    FMUModel(ME1/ME2).simulate(self, 
                        start_time=0.0,
                        final_time=1.0,
                        input=(),
                        algorithm='AssimuloAlg', 
                        options={}):

A list of the different solver specific parameters can be found in the :doc:`usage` section.
