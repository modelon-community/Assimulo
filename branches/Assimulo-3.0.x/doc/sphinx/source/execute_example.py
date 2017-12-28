from assimulo.solvers import *
from assimulo.problem import *
from assimulo.examples import *
from assimulo import examples
import os
import sys
import pylab as P

ex=sys.argv[1]
name=ex.split(os.sep)
script_name=name[-1]
plot_name=script_name.replace('.py','.png')  
name=os.path.join(os.getcwd(),plot_name)
reference_text='Assimulo Example:  {}'.format(script_name)

def savefig(filename = name):
    fig = P.gcf()
    fig.set_size_inches(*(1.0*fig.get_size_inches()))
    fig.set_facecolor([0.69,0.76,0.77])
    fig.text(0.05,0.03,reference_text,size='small')
    return fig.savefig(filename, facecolor=fig.get_facecolor())
P.show = savefig

examp = getattr(examples, script_name[:-3])
examp.run_example()
