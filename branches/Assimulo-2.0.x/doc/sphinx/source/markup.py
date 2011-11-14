from assimulo.solvers import *
from assimulo.problem import *

def mark_examples():
    import assimulo.examples as examples
    
    for ex in examples.__all__:
        file = open("EXAMPLE_"+ex+".rst",'w')
        
        file.write(ex + '.py\n')
        file.write('===================================\n\n')
        file.write('.. autofunction:: assimulo.examples.'+ex+'.run_example\n\n')
        file.write('.. note::\n\n')
        file.write('    Press [source] (to the right) to view the example.\n')
        file.close()
        
    

def mark_solvers():

    solvers = [(sundials.CVode, "ODE"), (sundials.IDA, "DAE"), (radau5.Radau5ODE, "ODE"), (radau5.Radau5DAE, "DAE"),
               (euler.ExplicitEuler, "ODE"), (runge_kutta.RungeKutta4, "ODE"), (runge_kutta.RungeKutta34, "ODE")]
    
    
    rhs = lambda t,y: [1.0]
    res = lambda t,y,yd: [1.0]
    
    exp_mod = Explicit_Problem(rhs, 0.0)
    imp_mod = Implicit_Problem(res, 0.0, 0.0)
    
    for solver in solvers:
        if solver[1] == "ODE":
            method = solver[0](exp_mod)
            str_ret = "t, y = "
        else:
            str_ret = "t, y, yd = "
            method = solver[0](imp_mod)
            
        options = method.get_options()
        supports = method.supports
        
        module_name = method.__class__.__module__.split(".")[-1]
        solver_name = method.__class__.__name__
        problem_name = method.problem.__class__.__name__
        
        file = open(solver[1]+"_"+solver_name+".rst",'w')
        
        file.write('\n')
        file.write(solver_name + '\n')
        file.write('=================================\n\n')
        file.write('Support\n----------------\n\n')
        file.write('- State events (root funtions) : '+str(supports["state_events"])+'\n')
        file.write('- Step events (completed step) : '+str(supports["step_events"])+'\n')
        file.write('- Time events : '+'True\n')
        file.write('\nUsage\n--------------\n\n')
        file.write('Import the solver together with the correct problem:: \n\n')
        file.write('    from assimulo.solvers.'+module_name+' import '+ solver_name+'\n')
        file.write('    from assimulo.problem import '+problem_name+'\n\n')
        file.write('Define the problem. \n\n')
        file.write('    :class:`'+problem_name+ ' <assimulo.problem.'+problem_name+'>`\n\n')
        file.write('Create and modify the solver parameters.\n\n')
        file.write('Parameters:\n\n')
        
        iter_options = options.keys()
        iter_options.sort()
        for opt in iter_options:

            str_name = '- :class:`' + opt + ' <assimulo.solvers.' + module_name + '.' +solver_name + '.'+ opt + '>`'
            
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
        file.write('    :class:`'+str_ret+solver_name+'.simulate(tfinal, ncp) <assimulo.solvers.'+module_name+'.'+solver_name+'.simulate>` \n\n')
        
        #file.write('Plot the solution.\n\n')
        #file.write('    :class:`'+solver_name+'.plot() <assimulo.solvers.'+module_name+'.'+solver_name+'.plot>`\n\n')
        
        file.write('Information:\n\n')
        file.write('- :class:`'+solver_name+'.get_options() <assimulo.solvers.'+module_name+'.'+solver_name+'.get_options>` Returns the current solver options.\n')
        file.write('- :class:`'+solver_name+'.get_supports() <assimulo.solvers.'+module_name+'.'+solver_name+'.get_supports>` Returns the functionality which the solver supports.\n')
        file.write('- :class:`'+solver_name+'.get_statistics() <assimulo.solvers.'+module_name+'.'+solver_name+'.get_statistics>` Returns the run-time statistics (if any).\n')
        file.write('- :class:`'+solver_name+'.print_statistics() <assimulo.solvers.'+module_name+'.'+solver_name+'.print_statistics>` Prints the run-time statistics for the problem.\n')
        
        file.close()

mark_examples()
mark_solvers()

