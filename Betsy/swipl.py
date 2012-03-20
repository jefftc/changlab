#Code to interact python and prolog using pexpect module. To run the
# method in this scripts, pexpect module and swi prolog should be installed.
# connect() method is to create a prolog handle.
# send_query() method is to send a prolog command and return a list of facts.
# source_code() method is to load a source pl code.
# disconnect() method is to close the prolog and disconnect the handle with
#            proglog.

import os
import pexpect
import string
import sys
os.environ["TERM"] = "dumb"
os.environ["TERMCAP"] = ""
waitingtime=16
def connect(plcmd=None):
    """Create a pexpect handle."""
    if plcmd == None:
        if sys.platform == "linux2":
            plcmd = "swipl"
        elif sys.platform == 'darwin':
            plcmd = "/opt/local/bin/swipl"
    p=pexpect.spawn(plcmd)
    return p

def send_query(handle,query):
    """Send query to prolog and return a list of facts."""
    if not query.endswith('.'):
        query = query + '.'
    handle.send(query+'\n')
    result = ''
    
    i = handle.expect(['true','false.','ERROR','='],timeout=waitingtime)
    #i=0, end of one solution,
    #i=1, no solution find,
    #i=2, Error when trying to find solution,
    #i=3, one solution end with variable assignment  
    if i == 1:
        result='' #if false, result is empty
        
    if i == 2:
        handle.expect('[#\?]',timeout=waitingtime)
        raise ValueError('ERROR'+handle.before) #if Error, raise error message
    if i == 3: # if variable assignment, find more solutions
        result = result + handle.before + '='
        handle.send(';')
        i = handle.expect(['true','false','ERROR','[#\?]',';'],
                                  timeout=waitingtime)
        if i not in [0,4]:   #if there is only one variable,will find 'false', need to get the result
            result = result + handle. before
    
    while i == 0 or i == 4:
        if i == 0:
            result = result + handle.before
            try: #to check if the result is "true.", if yes, no more solution,exit
                i = handle.expect(['[#\?]'],timeout=waitingtime)
                break   
            except pexpect.TIMEOUT:  # if TimeoutError, there are more solutions
                handle.send(';')
                i=handle.expect(['true','false','ERROR','[#\?]',';'],
                                timeout=waitingtime)
        if i == 4: #if found ';', to see if this is the last solution,if not, find more solutions
            result = result + handle.before + ';'
            try:
                i = handle.expect(['true','false','ERROR','[#\?]',';'],
                                  timeout=waitingtime)
            except pexpect.TIMEOUT:
                handle.send(';')
                i=handle.expect(['true','false','ERROR','[#\?]',';'],
                               timeout=waitingtime)
    result = result.split(';')
    newresult = [x for x in result if x.strip()]#get rid of the last empty line
    
    if newresult == []:
        return newresult 
    else:
        text = newresult[0].split('\r\n')  #clean the header
        newtext = [x.strip() for x in text if x.strip() and x!='.']
        newresult[0] = '\r\n'.join(newtext[1:])
        return newresult
    
        
def source_code(handle,sourcecode):
    """Load the pl source file."""
    if sourcecode.endswith('.pl'):
        sourcecode = sourcecode[:-3]
    assert os.path.exists(os.path.join(os.getcwd(),sourcecode+'.pl')
                          ),"Cannot find the source code"
    handle.send('['+sourcecode+'].'+'\n')
    loaded = handle.expect(['ERROR','true'],timeout=2)
    if loaded == 0:
         handle.expect('[#\?]',timeout=1)
         raise ValueError('ERROR'+handle.before)
    
def disconnect(handle):
    """Exit the prolog application and disconnect the father
       application to the child application.
       """
    handle.send('halt.'+'\n')
    handle.close()
    
def parse_solution(solution_str):
    text=''
    variables = dict()
    lines = solution_str.split('\r\n')
    newlines = lines[:]
    variable_line = []
    for line in lines:  #find the text line and variable line
        if '=' in line:
            newlines.remove(line)
            variable_line.append(line)
    text = '\r\n'.join(newlines)
    #find the '=' index in the variables line
    for i in variable_line:
        index_list = []
        start_list = [0]
        end_list = []
        start = -1
        while start < len(i)-1:
            try:
               start = i.index('=',start+1)
               index_list.append(start)
            except ValueError:
               break
        #if more than one in the same line, find the seperate point
        if len(index_list)>1:
            for j in range(len(index_list)-1):
                #find the last ',' between two '='
                seperate_point=i.rfind(',',index_list[j],index_list[j+1])
                end_list.append(seperate_point)
                start_list.append(seperate_point+1)
        end_list.append(len(i))
        for j in range(len(index_list)):
            variables[i[start_list[j]:index_list[j]].strip()] = i[index_list[j]+1:end_list[j]].strip()
    return text,variables
