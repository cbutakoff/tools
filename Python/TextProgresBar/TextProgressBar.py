# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 18:00:43 2017

@author: costa
needs python3
"""

"""
adapted from https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/34325723
Call in a loop to create terminal progress bar
@params:
    iteration   - Required  : current iteration (Int)
    total       - Required  : total iterations (Int)
    prefix      - Optional  : prefix string (Str)
    suffix      - Optional  : suffix string (Str)
    decimals    - Optional  : positive number of decimals in percent complete (Int)
    length      - Optional  : character length of bar (Int)
    fill        - Optional  : bar fill character (Str)
"""

class TextProgressBar:
    
    def __init__(self, total, prefix = '', suffix = '', decimals = 1, length = 50, fill = '█'):
        """
        @params:
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
        """
        self.prefix = prefix
        self.suffix = suffix
        self.decimals = decimals
        self.length = length
        self.fill = fill
        self.iteration = 0
        self.total = total
        self.filledLength = 0
        self.iteration = 0
        
    def UpdateProgress(self, iteration):
        """
        @params:
            iteration   - Required  : current iteration (Int)
        """
        self.iteration = iteration
        filledLength = int(self.length * iteration // self.total)
        if filledLength != self.filledLength:
            self.filledLength = filledLength
            self.__print()
        
    def __print(self):
        percent = ("{0:." + str(self.decimals) + "f}").format(100 * (self.iteration / float(self.total)))
        
        bar = self.fill * self.filledLength + '-' * (self.length - self.filledLength)
        print('\r%s |%s| %s%% %s' % (self.prefix, bar, percent, self.suffix), end = '\r')
        # Print New Line on Complete
        if self.iteration == self.total: 
            print()
        
        
        
        
def printProgressBar (iteration, total, ):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()
        