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
import inspect
from assimulo import testattr
from assimulo.special_systems import Mechanical_System
from assimulo.problem import Implicit_Problem
from assimulo.exception import *
from assimulo.lib import mexax
from scipy import *
from assimulo.tests.solvers import test_mexax_standalone
import unittest



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
            def mass_matrix(t,p):
                return eye(2)
            return Mechanical_System(n_p, forces, n_la,
                                         [1.,0.], [0.,0.],
                                         [0],
                                         [0.,0.], [0.,-g], GT=GT,
                                         constr3 = constr3,
                                         constr2 = constr2,
                                         constr1 = constr1,
                                         mass_matrix=mass_matrix)
        my_pend_sys=pendulum()
        self.index=index='oproj2'
        self.my_pend=my_pend_sys.generate_problem(index)
        self.my_pend.name='Index = {}'.format(index)
        # Build no the standalone solution which does not use Assimulo
        self.mexax_sta_test=test_mexax_standalone.Mexaxpendulum(2,gr=13.750371636040738)
    @testattr(stddist = True)
    def test_problem_fprob_mass(self):
        # compute standalone fprob
        qflag=9*[False]
        qflag[0]=True
        parameterlist={'nl':1,'ng':0,'nu':0,'t':0.,
                                  'p':array([0.,1.]),'v':array([0.,0.]),'u':array([0.]),
                                  'lam':array([0.]),'mass':empty((2,2)),'gp':empty((1,2)),
                                  'f':empty((2,)),'pdot':empty((2,)),'udot':empty((1,)),
                                  'g':empty((1,)),'gi':empty((1,)),'fl':empty((1,)),
                                  'qflag':qflag}
        self.mexax_sta_test.fprob(**parameterlist)
        self.my_pend.fprob(**parameterlist)
        assert((self.mexax_sta_test.fprob(**parameterlist)[0]==self.my_pend.fprob(**parameterlist)[0]).all())
        qflag[0]=False
        qflag[3]=True
        assert((self.mexax_sta_test.fprob(**parameterlist)[3]==self.my_pend.fprob(**parameterlist)[3]).all())
        qflag[3]=False
        assertRaises()
        qflag[2]=True
        

