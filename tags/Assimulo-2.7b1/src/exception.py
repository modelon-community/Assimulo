#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# Copyright (C) 2010 Modelon AB
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


class AssimuloException(Exception):
    pass 

class Algebraic_Exception(AssimuloException):
    pass
    
class AssimuloRecoverableError(AssimuloException):
    pass

class TerminateSimulation(AssimuloException):
    pass
    
class TimeLimitExceeded(AssimuloException):
    pass

class DiscardValue(AssimuloException):
    pass

class Explicit_ODE_Exception(AssimuloException):
    pass
    
class ODE_Exception(AssimuloException):
    pass

class Implicit_ODE_Exception(AssimuloException):
    pass    
    
class Rodas_Exception(AssimuloException):
    pass
    
class Dopri5_Exception(AssimuloException):
    pass
    
class GLIMDA_Exception(AssimuloException):
    pass

class ODEPACK_Exception(AssimuloException):
    pass

class DASP3_Exception(AssimuloException):
    pass
class RKStarter_Exception(AssimuloException):
    pass
