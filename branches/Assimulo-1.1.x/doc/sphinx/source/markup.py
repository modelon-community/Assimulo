
import assimulo.explicit_ode as Exp
import assimulo.implicit_ode as Imp

base_classes = ['ODE','Explicit_ODE','Implicit_ODE','Sundials', 'Radau_Common']

ignore = ['is_disc', 'h']
sundials = ['IDA', 'CVode']
make_consistency = ['IDA']
radau = ['Radau5']

methods = {'IDA':['make_consistency', 'interpolate', 'interpolate_sensitivity'] , 'CVode':['interpolate'],
           'Radau5':['interpolate']}



def copydoc(solver, type):
    
    classes = [solver]
    classes_names = [solver.__name__]
    bas = solver.__bases__
    
    while True:
        for i in bas:
            classes.append(i)
            classes_names.append(i.__name__)
        bas = bas[0].__bases__
        
        if not bas:
            break;
    
    name = solver.__name__
    sname = solver.__name__
    name = 'solver' + '_' + type + '_' + name + '.rst'
    file = open(name,'w')
    list_namnes = []
    
    if 'Explicit_ODE' in classes_names:
        class_inher = 'Explicit_ODE'
        module_base = 'explicit_ode'
    else:
        class_inher = 'Implicit_ODE'
        module_base = 'implicit_ode'
    
    file.write('\n')
    file.write(solver.__name__ + '\n')
    file.write('=================================\n\n')
    file.write('Usage\n--------------\n\n')
    file.write('Import the solver together with correct problem:: \n\n')
    file.write('    from assimulo.'+module_base+' import '+ solver.__name__+'\n')
    file.write('    from assimulo.problem import '+class_inher[:8]+'_Problem\n\n')
    file.write('Define the problem. \n\n')
    file.write('    :class:`'+class_inher[:8]+'_Problem <assimulo.problem.'+class_inher[:8]+'_Problem>`\n\n')
    file.write('Create and modify the solver parameters.\n\n')
    file.write('Parameters:\n\n')
    
    for solver in classes:
    
        for (name, values) in solver.__dict__.items():
            if not values.__class__.__name__.find('property'):
                if not classes_names[0] in sundials and name in ignore:
                    pass
                else:
                    if name not in list_namnes:
                        str_name = '- :class:`' + name + ' <assimulo.' + module_base + '.' +classes[0].__name__ + '.'+ name + '>`'
                        list_namnes.append(name)
                        str_doc  = ' ' + values.__doc__[0:values.__doc__.find('.')] + '.\n'
                        str_doc = str_doc.split()
                        str_doc_join = ' '
                        t = str_doc_join.join(str_doc)
                        str = str_name + ' ' +t + '\n'
                        file.write(str)
    
    try:
        lst =  methods[sname]
        file.write('\nMethods:\n\n')
        for met in lst:
            str_name = '- :class:`' + met + ' <assimulo.' + module_base + '.' +classes[0].__name__ + '.'+ met + '>`'
            #str_doc  = ' ' + values.__doc__[0:values.__doc__.find('.')] + '.\n'
            #str_doc = str_doc.split()
            #str_doc_join = ' '
            #t = str_doc_join.join(str_doc)
            #str = str_name + ' ' +t + '\n'
            file.write(str_name + '\n')
            
    except KeyError:
        pass
    
    file.write('\nSimulate the problem.\n\n')
    file.write('    :class:`'+classes_names[0]+'.simulate(tfinal, ncp) <assimulo.'+module_base+'.'+classes_names[0]+'.simulate>` \n\n')
    file.write('Plot the solution.\n\n')
    file.write('    :class:`'+classes_names[0]+'.plot() <assimulo.'+module_base+'.'+classes_names[0]+'.plot>`\n\n')
    
    file.write('Information.\n\n')
    file.write('- :class:`'+classes_names[0]+'.print_statistics() <assimulo.'+module_base+'.'+classes_names[0]+'.print_statistics>` Prints the run-time statistics for the problem.\n')
    
    if classes_names[0] in sundials:
        
        file.write('- :class:`'+classes_names[0]+'.print_event_info() <assimulo.'+module_base+'.'+classes_names[0]+'.print_event_info>` Prints the event information.\n')
    else:
        file.write('\n.. note::\n\n')
        file.write('    Only IDA and CVode supports discontinuous systems.')
    file.close()


def run():
    
    exp = 'exp'
    imp = 'imp'
    solvers = [(Exp.CVode, exp), (Imp.IDA, imp), (Exp.Explicit_Euler, exp), (Exp.RungeKutta4 ,exp)\
                , (Exp.RungeKutta34, exp), (Exp.Radau5, exp), (Imp.Radau5, imp)]
    for i in solvers:
        copydoc(i[0], i[1])
    
if __name__=='__main__':
    run()

