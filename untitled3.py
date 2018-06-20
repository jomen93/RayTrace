#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 14:06:31 2018

@author: ashcat
"""

import accretionStructures.simpleAccDisk as st
rDataArray = [10,4,3,2,1,1,2,3,4,5,10, 9, 2, 6, 3]

disk = st.SimpleDisk(rDataArray)

diskImage = disk.Ch()

print(diskImage)