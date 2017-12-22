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

import sys

#Import all the solvers from the different modules
try:
    from .euler import ExplicitEuler
    from .euler import ImplicitEuler
except ImportError as ie:
    sys.stderr.write("Could not find " + str(ie) + "\n")
try:
    from .radau5 import Radau5ODE
    from .radau5 import Radau5DAE
    from .radau5 import _Radau5ODE
    from .radau5 import _Radau5DAE 
except ImportError as ie:
    sys.stderr.write("Could not find " + str(ie) + "\n")
try:
    from .sundials import IDA
    from .sundials import CVode
except ImportError as ie:
    sys.stderr.write("Could not find " + str(ie) + "\n")
try:
    from .kinsol import KINSOL
except ImportError as ie:
    sys.stderr.write("Could not find " + str(ie) + "\n")
try:
    from .runge_kutta import RungeKutta34
    from .runge_kutta import RungeKutta4
    from .runge_kutta import Dopri5
except ImportError as ie:
    sys.stderr.write("Could not find " + str(ie) + "\n")
try:
    from .rosenbrock import RodasODE
except ImportError as ie:
    sys.stderr.write("Could not find " + str(ie) + "\n")
try:
    from .odassl import ODASSL
except ImportError as ie:
    sys.stderr.write("Could not find " + str(ie) + "\n")
try:
    from .odepack import LSODAR
except ImportError as ie:
    sys.stderr.write("Could not find " + str(ie) + "\n")
try:
    from .radar5 import Radar5ODE
except ImportError as ie:
    sys.stderr.write("Could not find " + str(ie) + "\n")
try:
    from .dasp3 import DASP3ODE
except ImportError as ie:
    sys.stderr.write("Could not find " + str(ie) + "\n")
try:
    from .glimda import GLIMDA
except ImportError as ie:
    sys.stderr.write("Could not find " + str(ie) + "\n")
