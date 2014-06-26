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

import nose
from assimulo import testattr
from assimulo.solvers.sundials import *
from assimulo.special_systems import Mechanical_System
from assimulo.problem import Implicit_Problem
from assimulo.exception import *

class Test_MEXAX:
    def setUp(self):
        """
        This sets up the pendulum test case.
        """
        def pendulum():
            g=13.7503671
            n_p=2
            n_la=1
            def forces(t, p, v):
                return N.array([0.,-g])
            def GT(p):
                return N.array([p[0],p[1]]).reshape((2,1))
            def constr3(t,y):
                p=y[0:2]
                return N.array([p[0]**2+p[1]**2-1.])
            def constr2(t,y):
                p,v=y[0:2],y[2:4]
                return N.array([p[0]*v[0]+p[1]*v[1]])
            def constr1(t,y):
                p,v,la=y[0:2],y[2:4],y[4:5]
                return N.array([v[0]**2+v[1]**2 - la[0] * (p[0]**2 + p[1]**2) - p[1] * g])
            return ass.Mechanical_System(n_p, forces, n_la,
                                         [1.,0.], [0.,0.],
                                         [0],
                                         [0.,0.], [0.,-g], GT=GT,
                                         constr3 = constr3,
                                         constr2 = constr2,
                                         constr1 = constr1)
        my_pend_sys=pendulum()
        index='oproj2'
        my_pend=my_pend_sys.generate_problem(index)
        my_pend.name='Index = {}'.format(index)