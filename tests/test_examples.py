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
from assimulo.exception import AssimuloException
import assimulo.examples as examples

class Test_Examples:
    @testattr(stddist = True)
    def test_cvode_with_jac_sparse(self):
        try:
            examples.cvode_with_jac_sparse.run_example(with_plots=False)
        except AssimuloException:
            pass #Handle the case when SuperLU is not installed
    
    @testattr(stddist = True)
    def test_ida_with_user_defined_handle_result(self):
        examples.ida_with_user_defined_handle_result.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_radau5dae_time_events_f(self):
        examples.radau5dae_time_events.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_kinsol_basic(self):
        examples.kinsol_basic.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_kinsol_with_jac(self):
        examples.kinsol_with_jac.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_kinsol_ors(self):
        examples.kinsol_ors.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_with_preconditioning(self):
        examples.cvode_with_preconditioning.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_dasp3_basic(self):
        print("Currently not running test_dasp3_basic. Numerically unstable problem.")
        #examples.dasp3_basic.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_gyro(self):
        examples.cvode_gyro.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_basic(self):
        examples.cvode_basic.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_cvode_with_disc(self):
        examples.cvode_with_disc.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_with_initial_sensitivity(self):
        examples.cvode_with_initial_sensitivity.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_with_jac(self):
        examples.cvode_with_jac.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_with_jac_spgmr(self):
        examples.cvode_with_jac_spgmr.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_ida_with_jac_spgmr(self):
        examples.ida_with_jac_spgmr.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_with_parameters(self):
        examples.cvode_with_parameters.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_cvode_with_parameters_fcn(self):
        examples.cvode_with_parameters_fcn.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_cvode_with_parameters_modified(self):
        examples.cvode_with_parameters_modified.run_example(with_plots=False)

    @testattr(stddist = True)
    def test_euler_basic(self):
        examples.euler_basic.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_euler_with_disc(self):
        examples.euler_with_disc.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_rungekutta4_basic(self):
        examples.rungekutta4_basic.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_rungekutta34_basic(self):
        examples.rungekutta34_basic.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_rungekutta34_with_disc(self):
        examples.rungekutta34_with_disc.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_ida_with_disc(self):
        examples.ida_with_disc.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_ida_with_initial_sensitivity(self):
        examples.ida_with_initial_sensitivity.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_ida_with_jac(self):
        examples.ida_with_jac.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_ida_with_parameters(self):
        examples.ida_with_parameters.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_radau5ode_vanderpol_c(self):
        examples.radau5ode_vanderpol.run_example(with_plots=False)
    
    @testattr(stddist = True)
    def test_radau5ode_vanderpol_f(self):
        examples.radau5ode_vanderpol.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_radau5ode_with_disc_c(self):
        examples.radau5ode_with_disc.run_example(with_plots=False)

    @testattr(stddist = True)
    def test_radau5ode_with_disc_f(self):
        examples.radau5ode_with_disc.run_example(with_plots=False)

    @testattr(stddist = True)
    def test_radau5ode_with_disc_sparse(self):
        examples.radau5ode_with_disc_sparse.run_example(with_plots=False)

    @testattr(stddist = True)
    def test_radau5ode_with_jac_sparse_c(self):
        examples.radau5ode_with_jac_sparse.run_example(with_plots=False)

    @testattr(stddist = True)
    def test_radau5dae_vanderpol_f(self):
        examples.radau5dae_vanderpol.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_dopri5_basic(self):
        examples.dopri5_basic.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_dopri5_with_disc(self):
        examples.dopri5_with_disc.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_rodasode_vanderpol(self):
        examples.rodasode_vanderpol.run_example(with_plots=False)
                
    @testattr(stddist = True)
    def test_mech_system_pendulum1(self):
        """
        This tests the class Mechanical_system together with ind1 and ida
        """
        examples.mech_system_pendulum.run_example('ind1',with_plots=False,with_test=True)
 
    @testattr(stddist = True)
    def test_mech_system_pendulum2(self):
        """
        This tests the class Mechanical_system together with ind2 and ida
        """
        examples.mech_system_pendulum.run_example('ind2',with_plots=False,with_test=True)
        
    @testattr(stddist = True)
    def test_mech_system_pendulum3(self):
        """
        This tests the class Mechanical_system together with ind3 and ida
        """
        examples.mech_system_pendulum.run_example('ind3',with_plots=False,with_test=True)
        
    @testattr(stddist = True)
    def test_mech_system_pendulum_ggl2(self):
        """
        This tests the class Mechanical_system together with ggl2 and ida
        """
        examples.mech_system_pendulum.run_example('ggl2',with_plots=False,with_test=True)
    
    @testattr(stddist = True)
    def test_mech_system_pendulum_ovstab2(self):
        """
        This tests the class Mechanical_system together with ovstab2 and odassl
        """
        examples.mech_system_pendulum.run_example('ovstab2',with_plots=False,with_test=True)
            
    @testattr(stddist = True)
    def test_mech_system_pendulum_ovstab1(self):
        """
        This tests the class Mechanical_system together with ovstab1 and odassl
        """
        examples.mech_system_pendulum.run_example('ovstab1',with_plots=False,with_test=True)
    
    @testattr(stddist = True)
    def test_lsodar_vanderpol(self):
        examples.lsodar_vanderpol.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_lsodar_with_disc(self):
        examples.lsodar_with_disc.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_euler_vanderpol(self):
        examples.euler_vanderpol.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_cvode_basic_backward(self):
        examples.cvode_basic_backward.run_example(with_plots=False)
        
    @testattr(stddist = True)
    def test_ida_basic_backward(self):
        examples.ida_basic_backward.run_example(with_plots=False)
