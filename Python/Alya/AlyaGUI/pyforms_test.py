import pyforms
from   pyforms          import BaseWidget
from   pyforms.controls import ControlCombo, ControlTextArea, ControlDir
from   pyforms.controls import ControlButton, ControlList, ControlText

import datetime
import os
import sys
import re

class AlyaFileIO():
    def __init__(self):
        self._casepath = './'
        self._is_ssh = False
        self._log_function = None  

    def SetLogControl(self, log_function):
        self._log_function = log_function

    def SetCasePath(self, path):
        self._casepath = path
        
    def GetCasePath(self):
        return self._casepath
    
    #def __OpenSSHSession(self):
    
        
    #def __CloseSSHSession(self):
    #lists all dat files
    def ListDats(self): 
        filelist = [];        
        
        if not self._is_ssh:            
            with os.scandir(self._casepath) as listOfEntries:  
                for entry in listOfEntries:
                    if entry.is_file():
                        if entry.name[-4:].lower()=='.dat':
                            filelist = filelist+[entry.name]
                       
        
        return filelist
        
    
    def ListCases(self):
        caselist = [];        
        
        if not self._is_ssh:            
            with os.scandir(self._casepath) as listOfEntries:  
                for entry in listOfEntries:
                    if entry.is_file():
                        if entry.name[-4:].lower()=='.dat':
                            if len(entry.name.split('.'))==2:
                                caselist = caselist+[entry.name.split('.')[0]]
                       
        
        return caselist
        
    def ReadFullFile(self, filename):
        contents = []
        
        #print(f'Reading {os.path.join(self._casepath,filename)}')
        if self._log_function != None:
            self._log_function(f'Reading {os.path.join(self._casepath,filename)}')
                    
        
        if not self._is_ssh:            
            try:
                with open(os.path.join(self._casepath,filename), 'r') as f:
                    contents = f.readlines()
            except:
                error_message = f"Unexpected error: {sys.exc_info()[0]}"
                if self._log_function == None:
                    print(error_message)
                else:
                    self._log_function(error_message)
                pass
                
        return contents
    
    
    def FindSectionNames(self, data_lines):
        section_names = []
        
        for line in data_lines:
            parts = re.split('\s+', line.strip())

            if parts[0].lower()[0:4] == 'end_':
                section_names += [parts[0].lower()[4:]]
                
        return section_names
        
    def ExtractImportantParameters(self, casename):

        params = ['alya', 'run_type', 'time_interval', 'number_of_steps']
        params_dom = ['nodal_points', 'elements','types_of_elements','boundaries','domain_integration_points']        
        params_ker = ['density','viscosity','steps']

        datfile = f'{casename}.dat'
        domdatfile = f'{casename}.dom.dat'
        kerdatfile = f'{casename}.ker.dat'
        nsidatfile = f'{casename}.nsi.dat'

        file_param_pairs = [(datfile, params), (domdatfile, params_dom), (kerdatfile, params_ker)]

        result_params = {};

        try:
            #read dat file
            for pair in file_param_pairs:
                with open(os.path.join(self._casepath, pair[0]), 'r') as f:
                    for line in f:
                        stripped_line = re.sub(r"\s+", "", line)
                        parts = re.split('[,:=]', stripped_line)
#                        print('parts:', parts)
                        
                        name = parts[0].lower()
                        value = None
                        
                        if len(parts)>1:
                            value = ", ".join(parts[1:])

#                            print('name:', name)
#                            print('value:', value)
        
                            if name in pair[1]:
                                result_params[name] = value

            #read boundary conditions
            reading_bdry_conditions = False
            boundary_conditions = []
            with open(os.path.join(self._casepath, nsidatfile), 'r') as f:
                for line in f:
                    stripped_line = re.sub(r"\s+", "", line)
                    parts = re.split('[,:=]', stripped_line)
                    
                    name = parts[0].lower()
                    if name=='boundary_conditions':
                        reading_bdry_conditions = True
                    elif name=='end_boundary_conditions':
                        reading_bdry_conditions = False
                    else:
                        if reading_bdry_conditions:
                            boundary_conditions += [line]
                    
            result_params['nsi_boundary_conditions'] = ''.join(boundary_conditions)
                
        except:
            error_message = f"Unexpected error reading {datfile}: {sys.exc_info()[0]}"
            if self._log_function == None:
                print(error_message)
            else:
                self._log_function(error_message)
            pass
        
        return result_params
        
    
    def FindAllIncludedFiles(self, datfile):
        #list included files in the specified dat file
        
        included_files = []        
        
        try:
            with open(os.path.join(self._casepath, datfile), 'r') as f:
                data = f.readlines()
                
            
            section_names = self.FindSectionNames(data)
            self._log_function(f'Found following sections: {section_names}')
                
            section_hierarchy = []
            for line in data:
                parts = re.split('\s+', line.strip())

                
                name = re.sub(r"^\s+|[\s,]+$", "", parts[0].lower())  #eliminate whitespaces and trailing commas

                if name in section_names:  #remove alsoo commas and other stuff
                    section_hierarchy += [name]
                elif name[0:4] == 'end_':
                    del section_hierarchy[-1]
                    
                
                
                if parts[0].lower()=='include':
                    included_files += [(parts[1].strip(), '/'.join(section_hierarchy))]
                
        except:
            error_message = f"Unexpected error reading {datfile}: {sys.exc_info()[0]}"
            if self._log_function == None:
                print(error_message)
            else:
                self._log_function(error_message)
            pass
    
        return included_files
#    def ParseDat(self, filename):       
#        try: 
#            with open(filename, 'r') as f:
#                lines = f.readlines();
            
        


#############################################################3
#
#
#
#       GUI
#
#
#
#
#
#############################################################        


class AlyaGUI(BaseWidget):


    
    def __init__(self):
        super(AlyaGUI,self).__init__('Alya Case Editor')

        #Layer for IO handling
        self._fileio = AlyaFileIO()

        self._log         = ControlTextArea('Log')
        self._fileio.SetLogControl(self.__addlog)



        #Main form fields
        self._casepath     = ControlDir('Case path', default='/home/costa/Downloads/1111')

        self._fileio.SetCasePath(self._casepath._value)


        self._casename    = ControlCombo('Case')

        
        
        
        self._casePathScan = ControlButton('Scan Path')
        self._casePathScan.value = self.___casePathScan_pressed
        
        self._caseReload = ControlButton('Load Case')
        self._caseReload.value = self.___caseReload_pressed
        
        self._casepath.changed_event = self.__casePathChanged
        
        #General tab. Here everythong most important to check
        self._dat_casename = ControlText('Name','')
        self._dat_run_type = ControlText('Run type','')
        self._dat_time_interval = ControlText('Time interval to run','')
        self._dat_number_of_steps = ControlText('Number of steps to run','')

        self._dom_nodal_points = ControlText('Number of points','')
        self._dom_elements = ControlText('Number of all elements','')
        self._dom_number_of_boundary_faces = ControlText('Number of boundary faces','')
        self._dom_mesh_numbers_verify_button = ControlButton('Verify')

        self._dom_types_of_elements = ControlText('Types of elements','')
        self._dom_domain_integration_points = ControlText('Number of integration points','')
        self._ker_density = ControlText('Density','')        
        self._ker_viscosity = ControlText('Viscosity','')        
        self._ker_steps_to_save = ControlText('Every how many steps to save','')        
        self._nsi_boundary_conditions = ControlTextArea('Nastin boundary conditions')
        general_tab = [('_dat_casename','_dat_run_type'),('_dat_time_interval','_dat_number_of_steps'),\
                       ('_dom_nodal_points','_dom_elements','_dom_number_of_boundary_faces','_dom_mesh_numbers_verify_button'),\
                       ('_dom_types_of_elements','_dom_domain_integration_points'),\
                       ('_ker_density','_ker_viscosity'),\
                       '_ker_steps_to_save','_nsi_boundary_conditions']
        
        #Raw casefile tab
        self._FileEditor = ControlTextArea('File editor')
        self._FileEditorSaveButton = ControlButton("Save file")
        self._DatPicker = ControlList('Dat files')
        self._DatPicker.horizontal_headers = ['Filename','Relative Path']
        self._DatPicker.item_selection_changed_event = self.__DAT_selection_changed
        self._DatPicker.readonly = True
        
        self._IncludedFilePicker = ControlList('Included files')
        self._IncludedFilePicker.horizontal_headers = ['Filename','Section']
        self._IncludedFilePicker.readonly = True
        self._IncludedFilePicker.item_selection_changed_event = self.__INCLUDE_selection_changed
        
        
        self.formset = [ ('_casepath','_casePathScan','||','_casename','_caseReload'),\
                        {'General': general_tab, \
                         'Raw Casefile':[('_FileEditor','||','_DatPicker','||','_IncludedFilePicker'),\
                         '_FileEditorSaveButton']},'=','_log']    

        self._log.autoscroll = True
        self._log.readonly = True
        
        
    def ___fill_general_tab(self, params):
        self._dat_casename.value = params['alya']
        self._dat_run_type.value = params['run_type']
        self._dat_time_interval.value = params['time_interval']
        self._dat_number_of_steps.value = params['number_of_steps']
        self._dom_nodal_points.value = params['nodal_points']
        self._dom_elements.value = params['elements']
        self._dom_number_of_boundary_faces.value = params['boundaries']
        self._dom_types_of_elements.value = params['types_of_elements']
        self._dom_domain_integration_points.value = params['domain_integration_points']
        self._ker_density.value = params['density']
        self._ker_viscosity.value = params['viscosity']
        self._ker_steps_to_save.value = params['steps']
        self._nsi_boundary_conditions.value = params['nsi_boundary_conditions']

        
    def ___caseReload_pressed(self):
        #
        # Reload the case
        #
        self.__addlog(f'Loading case {self._casename.value}')
        
        params = self._fileio.ExtractImportantParameters( self._casename.value )
        self.__addlog(f'Loaded parameters {params}')
        self.___fill_general_tab(params)


    def ___casePathScan_pressed(self):
        #
        # Scan path for the cases, fill the listbox of cases
        # Not fully implmented. Does not support more than one case per path yet
        #
        self._fileio.SetCasePath(self._casepath._value)
        self.__addlog(f'Scanning {self._casepath._value} for cases')
        
        cases = self._fileio.ListCases()
        self.__addlog(f'Found the following cases {cases}')

        self._casename.clear()
        if cases!=[]:
            for case in cases:
                self._casename.add_item(case)
                
            self._casename.current_index = 0
            self.___caseReload_pressed()
            
            dats = self._fileio.ListDats()
            self._DatPicker.clear()

            for dat in dats:
                self._DatPicker += [dat, '.'] 
       
    
    def __casePathChanged(self):
        self.__addlog(f'Case path changed to: {self._casepath._value}')
        self.___casePathScan_pressed()
        return True
    
    def __addlog(self, string):
        self._log.__add__(f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}: {string}')
        
        

    def __INCLUDE_selection_changed(self):
        #
        # Activated when user clicks on an INCLUDED file
        #
        selection = self._IncludedFilePicker.get_currentrow_value()        
        print(selection[0])
        data = self._fileio.ReadFullFile( selection[0] )

        self._FileEditor.value = ''
        if data!=[]:
            self._FileEditor.value = ''.join(data)

        
    def __DAT_selection_changed(self):
        #
        # Activated when user clicks on a dat file
        #
        selection = self._DatPicker.get_currentrow_value()   
        data = self._fileio.ReadFullFile( os.path.join(selection[1],selection[0]) )
        
        includes = self._fileio.FindAllIncludedFiles(selection[0])
        self.__addlog(f'Found {includes} included files')
        
        self._IncludedFilePicker.clear()
        if includes!=[]:
            for include in includes:
                self._IncludedFilePicker += [include[0], include[1]] 
            
        
        
        self._FileEditor.value = ''
        if data!=[]:
            self._FileEditor.value = ''.join(data)

#Execute the application
if __name__ == "__main__":   pyforms.start_app( AlyaGUI ) #, geometry=(100, 100, 2000, 2000) )
