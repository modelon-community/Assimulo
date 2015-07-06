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

def testattr(**kwargs):
    """Add attributes to a test function/method/class.
    
    This function is needed to be able to add
      @attr(slow = True)
    for functions.
    
    """
    def wrap(func):
        func.__dict__.update(kwargs)
        return func
    return wrap

try:
    import os
    curr_dir = os.path.dirname(os.path.abspath(__file__));
    _fpath=os.path.join(curr_dir,'version.txt')
    with open(_fpath, 'r') as f:
        __version__=f.readline().strip()
        __revision__=f.readline().strip()
except:
    __version__ = "unknown"
    __revision__ = "unknown"
