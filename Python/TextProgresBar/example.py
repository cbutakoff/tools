# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 18:16:00 2017

@author: costa
"""

from TextProgressBar import TextProgressBar 
from time import sleep

# A List of Items
items = list(range(0, 57))
l = len(items)

# Initial call to print 0% progress
pb = TextProgressBar(l, prefix = 'Progress:', suffix = 'Complete', length = 5)
for i, item in enumerate(items):
    # Do stuff...
    sleep(0.1)
    # Update Progress Bar
    pb.UpdateProgress(i + 1)


