import configparser
import re
import sys
import os

ProblemName = sys.argv[1]

ParameterFile = "parameters.dat"
TemplatePath = "./template"

config = configparser.ConfigParser(comment_prefixes='$')
config.read(ParameterFile)

subst_dict = {}

def apply_command(command, arguments):
    result = -1
    if command=='countlines': #arguemnnt = filename
        result = sum(1 for line in open(arguments.strip()))

    return result

for section in config.sections():
    #print('==',section)
    for key in config[section]:
        value = config[section][key]
        if '$' in value:
            value = value.split('$')[0]
        #print(f'{key} = {value}')

        subst_dict[key] = value.strip()
        if '![' in value:
            #if a command in ![]! is found - apply it to the arguments to the right 
            g = re.search('\!\[(.+)\]\!(.+)', value)
            #get 2 groups
            command = g.groups()[0]
            arguments = g.groups()[1]
            subst_dict[key] = apply_command(command, arguments)





print("Subsititution list:")
for key,value in subst_dict.items():
    print(f'{key} = "{value}"')


#get the list of all {ProblemName}.*.dat files
datfiles = [f for f in os.listdir(TemplatePath) if f.endswith('.tpl')]


for datfile in datfiles:
    print("Generating file:", datfile[:-4])

    with open(os.path.join(TemplatePath,datfile),'r') as f:
        data = f.read()

    for key,value in subst_dict.items():
        data = re.sub("(?i){%s}" % (key), f"{value}", data)

    with open(datfile[:-4],'w') as f:
        f.write(data)
