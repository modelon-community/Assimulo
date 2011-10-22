
class AssimuloException(Exception):
    pass 

class TerminateSimulation(AssimuloException):
    pass

class DiscardValue(AssimuloException):
    pass

class Explicit_ODE_Exception(AssimuloException):
    pass
