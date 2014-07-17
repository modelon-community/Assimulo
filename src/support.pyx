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
cimport numpy as N

from collections import OrderedDict

realtype = N.float

def set_type_shape_array(var, datatype=realtype):
    """
    Helper function to convert a scalar or list to a 1D-array
    """
    return  N.array(var, dtype = datatype).reshape(-1,)

class Statistics:
    def __init__(self):
        self.statistics = OrderedDict()
        self.statistics_msg = OrderedDict()
        
    def __setitem__(self, key, value):
        if self.statistics[key] == -1:
            self.statistics[key] = 0
        self.statistics[key] = value
        
    def __getitem__(self, key):
        if self.statistics[key] == -1:
            return 0
        else:
            return self.statistics[key]
        
    def add_key(self, key, msg):
        self.statistics[key] = -1
        self.statistics_msg[key] = msg
        
    def print_stats(self):
        max_len_msg = 0
        for k in list(self.statistics.keys()):
            if self.statistics[k] == -1:
                continue
            if max_len_msg < len(self.statistics_msg[k]):
                max_len_msg = len(self.statistics_msg[k])
        
        print("")
        for k in list(self.statistics.keys()):
            if self.statistics[k] == -1:
                continue
            print(" %s %s: %d")%(self.statistics_msg[k], " "*(max_len_msg-len(self.statistics_msg[k])+1) ,self.statistics[k])
        
    def reset(self):
        """
        Resets the statistics (to zero) that has been previously used
        """
        for k in list(self.statistics.keys()):
            if self.statistics[k] > -1:
                self.statistics[k] = 0
            
    def full_reset(self):
        """
        Resets all statistic to -1 (i.e. not used).
        """
        for k in list(self.statistics.keys()):
            self.statistics[k] = -1
        
    def keys(self):
        return self.statistics.keys()
