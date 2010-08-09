
from Assimulo.Explicit_ODE import *
from Assimulo.Implicit_ODE import *

base_classes = ['ODE','Explicit_ODE','Implicit_ODE','Sundials']

ignore = ['is_disc']
sundials = ['IDA', 'CVode']
make_consistency = ['IDA']

methods = {'IDA':['make_consistency']}


def copydoc(solver):
    

    
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
    name = 'solver' + '_' + name + '.rst'
    file = open(name,'w')
    list_namnes = []
    
    if 'Explicit_ODE' in classes_names:
        class_inher = 'Explicit_ODE'
    else:
        class_inher = 'Implicit_ODE'
    
    file.write('\n')
    file.write(solver.__name__ + '\n')
    file.write('=================================\n\n')
    file.write('Usage\n--------------\n\n')
    file.write('Import the solver together with correct problem:: \n\n')
    file.write('    from Assimulo.'+class_inher+' import '+ solver.__name__+'\n')
    file.write('    from Assimulo.Problem import '+class_inher[:8]+'_Problem\n\n')
    file.write('Define the problem. \n\n')
    file.write('    :class:`'+class_inher[:8]+'_Problem <Assimulo.Problem.'+class_inher[:8]+'_Problem>`\n\n')
    file.write('Create and modify the solver parameters.\n\n')
    file.write('Parameters:\n\n')
    
    for solver in classes:
    
        for (name, values) in solver.__dict__.items():
            if not values.__class__.__name__.find('property'):
                if not classes_names[0] in sundials and name in ignore:
                    pass
                else:
                    if name not in list_namnes:
                        str_name = '- :class:`' + name + ' <Assimulo.' + class_inher + '.' +classes[0].__name__ + '.'+ name + '>`'
                        list_namnes.append(name)
                        str_doc  = ' ' + values.__doc__[0:values.__doc__.find('.')] + '.\n'
                        str_doc = str_doc.split()
                        str_doc_join = ' '
                        t = str_doc_join.join(str_doc)
                        str = str_name + ' ' +t + '\n'
                        file.write(str)
                        
    file.write('\nSimulate the problem.\n\n')
    file.write('    :class:`'+classes_names[0]+'.simulate(tfinal, ncp) <Assimulo.'+class_inher+'.'+classes_names[0]+'.simulate>` \n\n')
    file.write('Plot the solution.\n\n')
    file.write('    :class:`'+classes_names[0]+'.plot() <Assimulo.'+class_inher+'.'+classes_names[0]+'.plot>`\n\n')
    
    if classes_names[0] in sundials:
        
        file.write('Information.\n\n')
        file.write('- :class:`'+classes_names[0]+'.print_statistics() <Assimulo.'+class_inher+'.'+classes_names[0]+'.print_statistics>` Prints the run-time statistics for the problem.\n')
        file.write('- :class:`'+classes_names[0]+'.print_event_info() <Assimulo.'+class_inher+'.'+classes_names[0]+'.print_event_info>` Prints the event information.\n')
    else:
        file.write('.. note::\n\n')
        file.write('    Only IDA and CVode supports discontinuous systems.')
    file.close()


def run():
    
    solvers = [CVode, IDA, Explicit_Euler, RungeKutta4, RungeKutta34]
    for i in solvers:
        copydoc(i)
    
    
if __name__=='__main__':
    
    run()

