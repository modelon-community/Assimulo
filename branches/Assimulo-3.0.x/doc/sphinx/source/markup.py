from assimulo.solvers import *
from assimulo.problem import *
import os

def mark_examples():
    import assimulo.examples as examples
    
    for ex in examples.__all__:
        file = open("EXAMPLE_"+ex+".rst",'w')
        
        file.write(ex + '.py\n')
        file.write('=============================================\n\n')
        file.write('.. autofunction:: assimulo.examples.'+ex+'.run_example\n\n')
        file.write('=============================================\n\n')
        file.write('.. program-output::   python '+os.path.join(os.getcwd(),'execute_example.py')+' '+os.path.join(os.getcwd(), examples.__path__[0]+os.sep+ex+'.py')+' \n\n')
        file.write('.. image:: '+os.sep+os.path.join(os.getcwd(),ex+'.png')+'\n\n')
        file.write('.. note::\n\n')
        file.write('    Press [source] (to the top right) to view the example.\n\n')
        file.close()
        
def mark_solvers():

    solvers = [(sundials.CVode, "ODE"), (sundials.IDA, "DAE"), (radau5.Radau5ODE, "ODE"), (radau5.Radau5DAE, "DAE"),
               (euler.ExplicitEuler, "ODE"), (runge_kutta.RungeKutta4, "ODE"), (runge_kutta.RungeKutta34, "ODE"),
               (runge_kutta.Dopri5, "ODE"), (rosenbrock.RodasODE, "ODE"), (odepack.LSODAR, "ODE"),(glimda.GLIMDA, "DAE"),
               (euler.ImplicitEuler, "ODE"), (dasp3.DASP3ODE, "ODE_SING"), (odassl.ODASSL,"DAE_OVER")]
    
    
    rhs = lambda t,y: [1.0]
    res = lambda t,y,yd: [1.0]
    dydt = lambda t,y,z: [1.0]
    dzdt = lambda t,y,z: [1.0]
    
    exp_mod = Explicit_Problem(rhs, 0.0)
    imp_mod = Implicit_Problem(res, 0.0, 0.0)
    sing_mod = SingPerturbed_Problem(dydt, dzdt, 0.0, 0.0)
    over_mod = Overdetermined_Problem(res,0.0,0.0)
    
    for solver in solvers:
        if solver[1] == "ODE":
            method = solver[0](exp_mod)
            str_ret = "t, y = "
        elif solver[1] == "DAE":
            str_ret = "t, y, yd = "
            method = solver[0](imp_mod)
        elif solver[1] == "ODE_SING":
            str_ret = "t, y = "
            method = solver[0](sing_mod)
        elif solver[1] == "DAE_OVER":
            str_ret = "t, y, yd = "
            method = solver[0](over_mod)
            
        options = method.get_options()
        supports = method.supports
        
        module_name = method.__class__.__module__.split(".")[-1]
        solver_name = method.__class__.__name__
        problem_name = method.problem.__class__.__name__
        
        file = open(solver[1]+"_"+solver_name+".rst",'w')
        
        file.write('\n')
        file.write(solver_name + '\n')
        file.write('=================================\n\n')
        file.write(solver[0].__doc__.replace("\n    ", "\n")) #REMOVES EXTRA INDENTATION
        file.write('\nSupport\n----------------\n\n')
        file.write('- State events (root funtions) : '+str(supports["state_events"])+'\n')
        file.write('- Step events (completed step) : '+str(supports["report_continuously"])+'\n')
        file.write('- Time events : '+'True\n')
        file.write('\nUsage\n--------------\n\n')
        file.write('Import the solver together with the correct problem:: \n\n')
        file.write('    from assimulo.solvers import '+ solver_name+'\n')
        file.write('    from assimulo.problem import '+problem_name+'\n\n')
        file.write('Define the problem, such as:: \n\n')
        
        if solver[1] == "ODE":
            file.write('    def rhs('+str_ret[:-3]+'): #Note that y are a 1-D numpy array.\n')
            file.write('        yd = -1.0\n')
            file.write('        return N.array([yd]) #Note that the return must be numpy array, NOT a scalar.\n\n')
            file.write('    y0 = [1.0]\n')
            file.write('    t0 = 1.0\n\n')
        elif solver[1] == "DAE":
            file.write('    def res('+str_ret[:-3]+'): #Note that y and yd are 1-D numpy arrays.\n')
            file.write('        res = yd[0]-1.0\n')
            file.write('        return N.array([res]) #Note that the return must be numpy array, NOT a scalar.\n\n')
            file.write('    y0  = [1.0]\n')
            file.write('    yd0 = [1.0]\n')
            file.write('    t0  = 1.0\n\n')
        elif solver[1] == "ODE_SING":
            file.write('    def rhs_slow(t,y,z): #Note that y and z are 1-D numpy arrays.\n')
            file.write('        return N.array([1.0]) #Note that the return must be numpy array, NOT a scalar.\n\n')
            file.write('    def rhs_fast(t,y,z): #Note that y and z are 1-D numpy arrays.\n')
            file.write('        return N.array([1.0]) #Note that the return must be numpy array, NOT a scalar.\n\n')
            file.write('    yy0 = [1.0]\n')
            file.write('    zz0 = [1.0]\n')
            file.write('    t0 = 1.0\n\n')
        elif solver[1] == "DAE_OVER":
            file.write('    def res('+str_ret[:-3]+'): #Note that y and yd are 1-D numpy arrays.\n')
            file.write('        res = [yd[0]-1.0, y[0]-1.0] \n')
            file.write('        return N.array([res]) #Note that the return must be numpy array, NOT a scalar.\n\n')
            file.write('    y0  = [1.0]\n')
            file.write('    yd0 = [1.0]\n')
            file.write('    t0  = 1.0\n\n')
    
        file.write('Create a problem instance::\n\n')
        if solver[1] == "ODE":
            file.write('    mod = '+problem_name+'(rhs, y0, t0)\n\n')
        elif solver[1] == "DAE":
            file.write('    mod = '+problem_name+'(res, y0, yd0, t0)\n\n')
        elif solver[1] == "ODE_SING":
            file.write('    mod = '+problem_name+'(rhs_slow, rhs_fast, yy0, zz0, t0)\n\n')
        elif solver[1] == "DAE_OVER":
            file.write('    mod = '+problem_name+'(res, y0, yd0, t0)\n\n')
        else:
            print "Unknown solver type"
        file.write('.. note::\n\n')
        file.write('    For complex problems, it is recommended to check the available :doc:`examples <examples>` and the documentation in the problem class, :class:`'+problem_name+ ' <assimulo.problem.'+problem_name+'>`. It is also recommended to define your problem as a subclass of :class:`'+problem_name+ ' <assimulo.problem.'+problem_name+'>`.\n\n')
        file.write('.. warning::\n\n')
        file.write('    When subclassing from a problem class, the function for calculating the right-hand-side (for ODEs) must be named *rhs* and in the case with a residual function (for DAEs) it must be named *res*.\n\n')
        file.write('Create a solver instance::\n\n')
        file.write('    sim = '+solver_name+'(mod)\n\n')
        file.write('Modify (optionally) the solver parameters.\n\n')
        file.write('    Parameters:\n\n')
        
        iter_options = options.keys()
        iter_options.sort()
        for opt in iter_options:

            str_name = '    - :class:`' + opt + ' <assimulo.solvers.' + module_name + '.' +solver_name + '.'+ opt + '>`'
            
            def find_doc(solv, opt):
                try:
                    str_doc = " ".join(solv.__dict__[opt].__doc__[:solv.__dict__[opt].__doc__.find(".")].split())
                except KeyError:
                    str_doc = " "
                    if len(solv.mro()) > 1:
                        str_doc = find_doc(solv.mro()[1], opt)
                        if str_doc == " " and len(solv.mro()) > 2:
                            str_doc = find_doc(solv.mro()[2], opt)
                return str_doc
            str_doc = find_doc(solver[0], opt)
            file.write(str_name + ' ' + str_doc + '.\n')
        
        if supports["interpolated_output"] or supports["interpolated_sensitivity_output"]:
            file.write('\nMethods:\n\n')
            file.write('- :class:`'+solver_name+'.interpolate <assimulo.solvers.'+module_name+'.'+solver_name+'.interpolate>`\n')
            
            if supports["interpolated_sensitivity_output"]:
                file.write('- :class:`'+solver_name+'.interpolate_sensitivity <assimulo.solvers.'+module_name+'.'+solver_name+'.interpolate_sensitivity>`\n')
            
            #SPECIAL FOR IDA
            if solver_name == "IDA":
                file.write('- :class:`'+solver_name+'.make_consistent <assimulo.solvers.'+module_name+'.'+solver_name+'.make_consistent>` Directs IDA to try to calculate consistent initial conditions.\n')
            
        file.write('\nSimulate the problem:\n\n')
        file.write('    :class:`'+str_ret+solver_name+'.simulate(tfinal, ncp, ncp_list) <assimulo.solvers.'+module_name+'.'+solver_name+'.simulate>` \n\n')
        
        #file.write('Plot the solution.\n\n')
        #file.write('    :class:`'+solver_name+'.plot() <assimulo.solvers.'+module_name+'.'+solver_name+'.plot>`\n\n')
        
        file.write('Information:\n\n')
        file.write('- :class:`'+solver_name+'.get_options() <assimulo.solvers.'+module_name+'.'+solver_name+'.get_options>` Returns the current solver options.\n')
        file.write('- :class:`'+solver_name+'.get_supports() <assimulo.solvers.'+module_name+'.'+solver_name+'.get_supports>` Returns the functionality which the solver supports.\n')
        file.write('- :class:`'+solver_name+'.get_statistics() <assimulo.solvers.'+module_name+'.'+solver_name+'.get_statistics>` Returns the run-time statistics (if any).\n')
        file.write('- :class:`'+solver_name+'.get_event_data() <assimulo.solvers.'+module_name+'.'+solver_name+'.get_event_data>` Returns the event information (if any).\n')
        file.write('- :class:`'+solver_name+'.print_event_data() <assimulo.solvers.'+module_name+'.'+solver_name+'.print_event_data>` Prints the event information (if any).\n')
        file.write('- :class:`'+solver_name+'.print_statistics() <assimulo.solvers.'+module_name+'.'+solver_name+'.print_statistics>` Prints the run-time statistics for the problem.\n')
        
        file.close()

mark_examples()
mark_solvers()

