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

import numpy as N
import pylab as P
import scipy as S
import scipy.linalg as LIN
import scipy.io as IO
import scipy.sparse as SPARSE
import scipy.sparse.linalg as LINSP
import nose
import os
from assimulo.solvers import KINSOL
from assimulo.problem import Algebraic_Problem
import warnings
import scipy.sparse


warnings.simplefilter("ignore", scipy.sparse.SparseEfficiencyWarning)

file_path = os.path.dirname(os.path.realpath(__file__))

def run_example(with_plots=True):
    r"""
    Example to demonstrate the use of the Sundials solver Kinsol with
    a user provided Jacobian and a preconditioner. The example is the 
    'Problem 4' taken from the book by Saad:
    Iterative Methods for Sparse Linear Systems.
    
    on return:
    
       - :dfn:`alg_mod`    problem instance
    
       - :dfn:`alg_solver`    solver instance
    
    """
    #Read the original matrix
    A_original = IO.mmread(os.path.join(file_path,"kinsol_ors_matrix.mtx"))

    #Scale the original matrix
    A = SPARSE.spdiags(1.0/A_original.diagonal(), 0, len(A_original.diagonal()), len(A_original.diagonal())) * A_original

    #Preconditioning by Symmetric Gauss Seidel
    if True:
        D = SPARSE.spdiags(A.diagonal(), 0, len(A_original.diagonal()), len(A_original.diagonal()))
        Dinv = SPARSE.spdiags(1.0/A.diagonal(), 0, len(A_original.diagonal()), len(A_original.diagonal()))
        E = -SPARSE.tril(A,k=-1)
        F = -SPARSE.triu(A,k=1)
        L = (D-E).dot(Dinv)
        U = D-F
        Prec = L.dot(U)
        
        solvePrec = LINSP.factorized(Prec)

    #Create the RHS
    b = A.dot(N.ones((A.shape[0],1)))
    
    #Define the res
    def res(x):
        return A.dot(x.reshape(len(x),1))-b
        
    #The Jacobian
    def jac(x):
        return A.todense()
    
    #The Jacobian*Vector
    def jacv(x,v):
        return A.dot(v.reshape(len(v),1))
    
    def prec_setup(u,f, uscale, fscale):
        pass
    
    def prec_solve(r):
        return solvePrec(r)
        
    y0 = S.rand(A.shape[0])
    
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
        solver.ftol = LIN.norm(res(solver.y0))*1e-9
        solver.max_iter = 300
        solver.verbosity = 10
        solver.globalization_strategy = "none"
        
    setup_param(alg_solver)
    setup_param(alg_solver_prec)
    
    #Solve orignal system
    y = alg_solver.solve()

    #Solve Preconditionined system
    y_prec = alg_solver_prec.solve()
    
    print("Error                 , in y: ", LIN.norm(y-N.ones(len(y))))
    print("Error (preconditioned), in y: ", LIN.norm(y_prec-N.ones(len(y_prec))))
    
    if with_plots:

        P.figure(4)
        P.semilogy(alg_solver.get_residual_norm_nonlinear_iterations(), label="Original")
        P.semilogy(alg_solver_prec.get_residual_norm_nonlinear_iterations(), label='Preconditioned')
        P.xlabel("Number of Iterations")
        P.ylabel("Residual Norm")
        P.title("Solution Progress")
        P.legend()
        P.grid()
        
        P.figure(5)
        P.plot(y, label="Original")
        P.plot(y_prec, label="Preconditioned")
        P.legend()
        P.grid()
        
        P.show()
    
    #Basic test
    for j in range(len(y)):
        nose.tools.assert_almost_equal(y[j], 1.0, 4)
        
    return [alg_mod, alg_mod_prec], [alg_solver, alg_solver_prec]

if __name__=='__main__':
    mod, solv = run_example()

