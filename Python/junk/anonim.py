import pydicom
import sys
import os
import shutil

path_in = sys.argv[1]
path_out = sys.argv[2]

def person_names_callback(dataset, data_element):
    if data_element.VR == "PN":
        data_element.value = "anonymous"

def anonymize_file(filenamein, filenameout):

    dataset = pydicom.dcmread(filenamein)


    data_elements = ['PatientID',
                     'PatientBirthDate']

    for de in data_elements:
        print(dataset.data_element(de))


    dataset.PatientID = "id"
    dataset.walk(person_names_callback)

    if 'OtherPatientIDs' in dataset:
        delattr(dataset, 'OtherPatientIDs')

    if 'OtherPatientIDsSequence' in dataset:
        del dataset.OtherPatientIDsSequence

    tag = 'PatientBirthDate'
    if tag in dataset:
        dataset.data_element(tag).value = '19000101'

    data_elements = ['PatientID',
                     'PatientBirthDate']
    for de in data_elements:
        print(dataset.data_element(de))

    dataset.save_as(filenameout)




#==================================================
for dirpath, dirnames, filenames in os.walk(path_in):
   for name in filenames:
      print(dirpath, ' == ',  name)
      infile = os.path.join(dirpath, name)
      outfile = os.path.join(path_out,dirpath, name)

      os.makedirs( os.path.join(path_out,dirpath), exist_ok=True )
      anonymize_file(infile, outfile)

