import datetime
import os
import sys
import re
import pprint
from collections import deque

#global variable to be added to the repeating keys
#e.g. repeating 'include' will become 'include%1', 'include%2',...
id_repeating_keys = 0;
    
class AlyaFileIO():
    def __init__(self):
        self._casepath = './'
        self._is_ssh = False

    def SetCasePath(self, path):
        self._casepath = path
        
    def GetCasePath(self):
        return self._casepath

    def GetFileSize(self, filename):
        size = -1
        
        if not self._is_ssh:
            try:            
                size = os.path.getsize( os.path.join(self._casepath,filename) )
            except:
                #print("Exception caught:", sys.exc_info()[0])
                #raise
                print('Error getting filesize of  ',os.path.join(self._casepath, filename))
        
        return size

    def ReadLines(self, datfile):
        lines = None
        
        if not self._is_ssh:
            try: 
                #with open(datfile, 'r') as f:
                #    lines = f.readlines()
                
                donotread = False;
                
                #do not read files that start with $NOPARSE
                with open(os.path.join(self._casepath, datfile)) as f:
                    line = re.sub('\s*','',f.readline().strip())
                    if '$NOPARSE' in line:
                        donotread = True
                     
                if not donotread:
                    lines = deque(open(os.path.join(self._casepath, datfile)))
                else:
                    print('File ',  os.path.join(self._casepath, datfile), ' is not parsed as per request')
            except:
                #print("Exception caught:", sys.exc_info()[0])
                #raise
                print('Error reading file ',os.path.join(self._casepath, datfile))
                    
        return lines
    
    def WriteLines(self, lines, datfile):
        if not self._is_ssh:
            try: 
                with open(os.path.join(self._casepath, datfile), 'w') as f:
                    f.writelines(lines)
            except:
                print("Exception caught:", sys.exc_info()[0])
                raise
        
        

class AlyaCaseParser():
    def __init__(self):
        self.__FileReader = AlyaFileIO()
        self.__IncludeFileSizeLimit = 100000 #do not read files bigger than 100k
        
        #module - file extension correspondence
        self.__ModuleFileCorrespondence = {'nastin_module':'nsi',  'temper_module':'tem', 'exmedi_module':'exm'}
        self.__CaseName = ''
        self.__Hierarchy = {}
        self.__NSpacesPerIndent = 3
        
    def SetCasePath(self, path):
        self.__FileReader.SetCasePath(path)
        
    def SetCaseName(self, casename):
        self.__CaseName = casename
        
    def GetHierarchy(self):
        return self.__Hierarchy
    
    def LoadCase(self):
        self.__Hierarchy = {}
        self.__Hierarchy['MAIN'] = self.__ParseDat(f'{self.__CaseName}.dat')
        self.__Hierarchy['KERNEL'] = self.__ParseDat(f'{self.__CaseName}.ker.dat')
        self.__Hierarchy['DOMAIN'] = self.__ParseDat(f'{self.__CaseName}.dom.dat')
        
        
        #find modules and read dat files
        for key, value in self.__Hierarchy['MAIN']['problem_data'].items():
            if key in self.__ModuleFileCorrespondence:
                self.__Hierarchy[key.upper()] = self.__ParseDat(f'{self.__CaseName}.{self.__ModuleFileCorrespondence[key]}.dat') 
    
        return self.__Hierarchy

    
    def __ParseDat(self, datfile):    
        #scan through the dat file and create a hierarchy of sections based n "end_" names
        
        hierarchy = {}
        
        lines = self.__FileReader.ReadLines(datfile)    
        hierarchy = self.__ParseSection(lines)
        
        return hierarchy
                    
    def __ReadIncludedFile(self, filename):
        lines = None

        if self.__FileReader.GetFileSize(filename)<self.__IncludeFileSizeLimit:
            lines = self.__FileReader.ReadLines(filename)    
        else:
            print('Included file ', filename, ' is too big to read')
        
        return lines
    
    
    def __CleanupLine(self, line):
        line = line.replace(':',' ')
        line = line.replace('=',' ')
        line = line.replace(',',' ')
        line = re.sub('\s*&\s*','&',line)
        return line
    
    #lines has to be a deque
    def __ParseSection(self, lines):
        hierarchy = {}
        global id_repeating_keys

        while lines: #parse backwards, while deque is not empty
            line = lines.pop()
            
            #find comments and split the lines into content and comments
            comment_index = line.find('$')
            
            comment = ''
            if comment_index>=0:  #skip comments at the end of the line
                comment = line[comment_index+1:].strip() 
                line = line[:comment_index].strip()
            
            if line!='': 
                line = self.__CleanupLine(line)
                parts = line.strip().split()

                #look for the sections from the end of file, sections are identified
                #by 'end_' statement at the beginning
                #if the section is found all the lines are buffered and resent to this function
                #otherwise process every line
                if 'end_' in parts[0].lower():  #if we have found a section, do a recursive call
                    section_name = parts[0].strip()[4:].lower()

                    #if this section name already exists, add a unique id
                    if section_name in hierarchy:
                        #if there are repeating keywords, use %X to add an id, for python                        
                        id_repeating_keys = id_repeating_keys+1                        
                        section_name_withid = section_name+f'%{id_repeating_keys}'
                    else:
                        section_name_withid = section_name


                    
                    line_subset = deque([]) #subset of lines stored here
                    
                    while lines: #copy out the whole section for recursive call

                        comment = ''
                        line = lines.pop()
                        #line = line.lower()
                        comment_index = line.find('$')
                        if comment_index>=0:  #skip comments at the end of the line
                            comment = line[comment_index+1:].strip() 
                            line = line[:comment_index].strip()


                        line = self.__CleanupLine(line)



                        parts1 = line.strip().lower().split()
                        if parts1==[]:
                            continue
                            
                        if parts1[0].strip()[0]=='$': #ignore comments
                            continue
                        
                        if parts1[0].strip()[0:len(section_name)] == section_name:

                            #keep stuff after section name, these are sectin options
                            #also remove spaces around '='
                            #if the section name had = and number glued to it, the right side of = 
                            #is added to the contents
                            pp = []
                            if len(parts1[0].strip()) > len(section_name):
                                pp = [parts1[0].strip()[len(section_name):].strip()]
                                
                            if len(parts1)>1:
                                pp += parts1[1:]
                            
                            options = pp
                            break
                        else:
                            #store lines in straight (not reversed) order
                            line_subset.appendleft(line)
                            options = []
                        
                    hierarchy[section_name_withid] = self.__ParseSection(line_subset) 
                    if(options!=[]):
                        if len(options)!=1:
	                        hierarchy[section_name_withid]['OPTIONS'] = options 
                        else:
	                        hierarchy[section_name_withid]['OPTIONS'] = options[0]
                    
                    if(comment!=''):
                        hierarchy[section_name_withid]['COMMENT'] = comment 
                        
                else: #otherwise at this level create a key with all the properties
                    name = parts[0].lower()

                    pp = []
                    if len(parts)>1:
                        pp += parts[1:]
                    
                    #included files are appended to the end of the deque at this level
                    #if they are too big, they stay as 'include':['filename'] pair
                    if name != 'include':      
                        if len(pp)!=1:                  
	                        hierarchy[ name ] = pp 
                        else:
	                        hierarchy[ name ] = pp[0] 
    
                        if(comment!=''):
                            hierarchy[name] += ['$COMMENT: '+str(comment)]
                    else:
                        included_lines = self.__ReadIncludedFile(parts[1].strip())
                        
                        if included_lines==None:
                            if name in hierarchy:
                                #if there are repeating keywords, use %X to add an id, for python
                                id_repeating_keys = id_repeating_keys+1                                
                                name = name+f'%{id_repeating_keys}'
                                
                            hierarchy[ name ] = pp
        
                            if(comment!=''):
                                hierarchy[name] += ['COMMENT: '+str(comment)]
                        else:
                            while included_lines:
                                lines.append(included_lines.popleft())
                            
                            


                    
        return hierarchy;                
        
    
    def SaveCase(self):
        self.__SaveHierarchy( self.__Hierarchy['MAIN'], f'{self.__CaseName}.dat')
        self.__SaveHierarchy( self.__Hierarchy['KERNEL'], f'{self.__CaseName}.ker.dat')
        self.__SaveHierarchy( self.__Hierarchy['DOMAIN'], f'{self.__CaseName}.dom.dat')
                
        #find modules and read dat files
        for key, value in self.__Hierarchy.items():
            if '_module' in key.lower():
                self.__SaveHierarchy( self.__Hierarchy[str(key)], \
                    f'{self.__CaseName}.{self.__ModuleFileCorrespondence[key.lower()]}.dat')
        
      
        
    def __HierarchyToLines(self, hierarchy, nidents):
        #idents - number of spaces to prepend
        lines = []
        
        nspaces = ' '*self.__NSpacesPerIndent*nidents;
        
        for key, value in hierarchy.items():
            if key!='OPTIONS':
                if type(hierarchy[key]) == dict:
                    options = ''
                    if 'OPTIONS' in hierarchy[key]:
                        options = hierarchy[key]['OPTIONS']

                    tt = key.split('%')

                    if type(options)==list:
	                    lines += [ nspaces + tt[0].upper()+' '+' '.join(options)+'\n' ]
                    else:
	                    lines += [ nspaces + tt[0].upper()+' '+options+'\n' ]

                    lines += self.__HierarchyToLines(hierarchy[key], nidents+1)
                    lines += [ nspaces + f'end_{tt[0]}'.upper() +'\n']
                else:
                    tt = key.split('%')
                    if type(value)==list:
                        lines += [ nspaces + tt[0].upper() +' ' + ' '.join(value) +'\n']
                    else:
                        lines += [ nspaces + tt[0].upper() +' ' + value +'\n']
        
        return lines
        
    def __SaveHierarchy(self, hierarchy, filename):
        print('save', filename)
        lines = self.__HierarchyToLines(hierarchy, 0)
        self.__FileReader.WriteLines(lines, filename)
        
        

                
pp = pprint.PrettyPrinter(indent=4)
            
parser = AlyaCaseParser()
parser.SetCasePath('/home/costa/Downloads/1111')
parser.SetCaseName('fluidda')
pp.pprint( parser.LoadCase() )
aa = parser.GetHierarchy()

parser.SetCasePath('/home/costa/Downloads/1111/save/')
parser.SaveCase()


import json
js = json.dumps(aa)
with open('jsonout.json','w') as f:
	f.write(js)

