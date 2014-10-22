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

__all__ = ["euler","radau5","sundials","runge_kutta","rosenbrock",
           "glimda","odepack","radar5","dasp3","odassl"]

#Import all the solvers from the different modules
from .euler import ExplicitEuler, ImplicitEuler
from .radau5 import Radau5ODE, Radau5DAE, _Radau5ODE, _Radau5DAE 
from .sundials import IDA, CVode
from .kinsol import KINSOL
from .runge_kutta import RungeKutta34, RungeKutta4, Dopri5
from .rosenbrock import RodasODE
from .odassl import ODASSL
from .odepack import LSODAR
from .radar5 import Radar5ODE
try:
    from .dasp3 import DASP3ODE
except ImportError as ie:
    print("Could not find {}".format(ie.args[0][16:]))
try:
    from .glimda import GLIMDA
except ImportError as ie:
    print("Could not find {}".format(ie.args[0][16:]))
