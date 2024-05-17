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

import os
import nose
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import scipy.io as spi
from assimulo.solvers import KINSOL
from assimulo.problem import Algebraic_Problem
import warnings

warnings.simplefilter("ignore", sps.SparseEfficiencyWarning)

file_path = os.path.dirname(os.path.realpath(__file__))

def run_example(with_plots=True):
    r"""
    Example to demonstrate the use of the Sundials solver Kinsol with
    a user provided Jacobian and a preconditioner. The example is the 
    'Problem 4' taken from the book by Saad:
    Iterative Methods for Sparse Linear Systems.
    """
    #Read the original matrix
    A_original = spi.mmread(os.path.join(file_path,"kinsol_ors_matrix.mtx"))

    #Scale the original matrix
    A = sps.spdiags(1.0/A_original.diagonal(), 0, len(A_original.diagonal()), len(A_original.diagonal())) * A_original

    #Preconditioning by Symmetric Gauss Seidel
    if True:
        D = sps.spdiags(A.diagonal(), 0, len(A_original.diagonal()), len(A_original.diagonal()))
        Dinv = sps.spdiags(1.0/A.diagonal(), 0, len(A_original.diagonal()), len(A_original.diagonal()))
        E = -sps.tril(A,k=-1)
        F = -sps.triu(A,k=1)
        L = (D-E).dot(Dinv)
        U = D-F
        Prec = L.dot(U)
        
        solvePrec = spsl.factorized(Prec)

    #Create the RHS
    b = A.dot(np.ones(A.shape[0]))
    
    #Define the res
    def res(x):
        return A.dot(x) - b
        
    #The Jacobian
    def jac(x):
        return A.todense()
    
    #The Jacobian*Vector
    def jacv(x,v):
        return A.dot(v)
    
    def prec_setup(u,f, uscale, fscale):
        pass
    
    def prec_solve(r):
        return solvePrec(r)
        
    y0 = np.random.rand(A.shape[0])
    
    #Define an Assimulo problem
    alg_mod = Algebraic_Problem(res, y0=y0, jac=jac, jacv=jacv, name = 'ORS Example')
    alg_mod_prec = Algebraic_Problem(res, y0=y0, jac=jac, jacv=jacv, prec_solve=prec_solve, prec_setup=prec_setup, name = 'ORS Example (Preconditioned)')

    #Define the KINSOL solver
    alg_solver = KINSOL(alg_mod)
    alg_solver_prec = KINSOL(alg_mod_prec)
    
    #Sets the parameters
    def setup_param(solver):
        solver.linear_solver = "spgmr"
        solver.max_dim_krylov_subspace = 10
        solver.ftol = np.linalg.norm(res(solver.y0))*1e-9
        solver.max_iter = 300
        solver.verbosity = 10
        solver.globalization_strategy = "none"
        
    setup_param(alg_solver)
    setup_param(alg_solver_prec)
    
    #Solve original system
    y = alg_solver.solve()

    #Solve Preconditioned system
    y_prec = alg_solver_prec.solve()
    
    print("Error                 , in y: ", np.linalg.norm(y-np.ones(len(y))))
    print("Error (preconditioned), in y: ", np.linalg.norm(y_prec-np.ones(len(y_prec))))
    
    if with_plots:
        import pylab as pl
        pl.figure(4)
        pl.semilogy(alg_solver.get_residual_norm_nonlinear_iterations(), label="Original")
        pl.semilogy(alg_solver_prec.get_residual_norm_nonlinear_iterations(), label='Preconditioned')
        pl.xlabel("Number of Iterations")
        pl.ylabel("Residual Norm")
        pl.title("Solution Progress")
        pl.legend()
        pl.grid()
        
        pl.figure(5)
        pl.plot(y, label="Original")
        pl.plot(y_prec, label="Preconditioned")
        pl.legend()
        pl.grid()
        
        pl.show()
    
    #Basic test
    for j in range(len(y)):
        nose.tools.assert_almost_equal(y[j], 1.0, 4)
        
    return [alg_mod, alg_mod_prec], [alg_solver, alg_solver_prec]

if __name__=='__main__':
    mod, solv = run_example()
