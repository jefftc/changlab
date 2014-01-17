# -*- coding: utf-8 -*-

upload_path = '/home/xchen/chencode/upload/'

def status():
    return dict(request=request,session=session,response=response)

@auth.requires_login()
def index():
    import os
    form_queue = FORM(INPUT(_type='submit',_value='Show Jobs',_name='status', _style ='height: 60px;width:120px;font-size: 20px;font-weight:bold'))
    result=dict(form_queue=form_queue)
    for label in labels:
        vars()['form_'+ label]=FORM(INPUT(_type='text',_value=label_showname[label],_name=label))
        result['form_'+ label]=vars()['form_'+ label]
    if form_queue.process(formname='status').accepted:
        redirect(URL('show_status'))
    return result
    

def protocol():
    x=get_protocol_info(request.args[0])
    protocol_parameters,protocol_parameters_dict,protocol_inputs,column_name,pretty_name = x
    fields = create_form(protocol_inputs,protocol_parameters)
    form=SQLFORM.factory(*fields)
    form_queue = FORM(INPUT(_type='submit',_value='show jobs status',_style ='height: 50px;width:200px;font-size: 18px;font-weight:bold'))
    if form.process().accepted:
        form_value = form.vars
        input_file,command_parameters = create_command(form_value)
        session.parameters = command_parameters
        session.input_file = input_file
        session.protocol= request.args[0]
        redirect(URL('submit_job'))
    if form_queue.process().accepted:
        redirect(URL('show_status'))
    return dict(form=form,form_queue=form_queue,protocol_inputs=protocol_inputs,
                protocol_parameters=protocol_parameters_dict,
                column_name=column_name,pretty_name=pretty_name)

def show_status():
    import subprocess
    all_cmd = ['python','/home/xchen/chencode/queue/queue.py','--list']    
    process = subprocess.Popen(all_cmd, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    process.wait()
    log_information = process.stdout.read().split('\n')[1:]
    log_information = [i for i in log_information if i]
    for log_line in log_information:
        lines = log_line.split('\t')
        db.status.update_or_insert(db.status.Jobnumber==lines[0],Jobnumber=lines[0],
                Jobname=lines[1], Status=lines[2], Started=lines[3],
                Ended=lines[4],Submitted=lines[5],Command=lines[6],User=lines[7])
    form=db().select(db.status.ALL)
    return dict(form=form)

    

def show_job():
    "show a job page"
    session.output=[]
    show_status()
    this_job = db.status(request.args(0)) or redirect(URL('index'))
    db.status.id.default=this_job.id
    text = ''
    newtext=''
    if  this_job.Status in ['running','completed']:
        hash_string = hash_command(this_job.Submitted,this_job.Command)
        pro_file = os.path.join('/home/xchen/chencode/queue/processing_info/', hash_string+'.txt')
        f=file(pro_file,'r')
        text = f.read().split('\n')
        f.close()
        newtext=text[:]
        for i in range(len(text)):
            if text[i].startswith('Report:'):
                text[i]=text[i][7:]
            if os.path.exists(text[i].strip()):
                newtext[i]=os.path.join(text[i].strip().split('/')[-2:])
            else:
                newtext[i]=text[i]
    info = 'Job is ' + this_job.Status
    return dict(job=this_job,info=info,text=text,newtext=newtext,session=session)

    

def submit_job():
    import subprocess
    import random
    import string
    jobname = ''.join(random.choice(string.ascii_uppercase+string.digits)
                          for x in range(6))
    cmd = ['python', '/home/xchen/chencode/Betsy/scripts/run_protocol.py', 
           '--protocol', session.protocol,'--user',str(auth.user.first_name+ auth.user.last_name),'--job_name',jobname] + session.input_file  + session.parameters
    all_cmd = ['python','/home/xchen/chencode/queue/queue.py' ,'--user',str(auth.user.first_name+ auth.user.last_name),'-j', jobname, ' ' .join(cmd)]
    process = subprocess.Popen(all_cmd, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)    
    process.wait()                      
    redirect(URL('show_status'))
    return dict()

    
def result_folder():
    args = '/'.join(request.args)
    path = os.path.join('/home/xchen/chencode/examples/sample_analysis',args)
    return response.stream(path)
    #return open(path,'rb').read()

def delete_job():
    jobname = request.args(0)
    cmd = ['python','/home/xchen/chencode/queue/queue.py','--kill',jobname]
    process = subprocess.Popen(cmd, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)    
    process.wait()   
    delete_job_and_update(jobname)
    redirect(URL('show_status'))
    
def user():
    """
    exposes:
    http://..../[app]/default/user/login
    http://..../[app]/default/user/logout
    http://..../[app]/default/user/register
    http://..../[app]/default/user/profile
    http://..../[app]/default/user/retrieve_password
    http://..../[app]/default/user/change_password
    use @auth.requires_login()
        @auth.requires_membership('group name')
        @auth.requires_permission('read','table name',record_id)
    to decorate functions that need access control
    """
    return dict(form=auth())

def download():
    """
    allows downloading of uploaded file
    http://..../[app]/default/download/[filename]
    """
    
    return response.download(request, db)

   
def call():
    """
    exposes services. for example:
    http://..../[app]/default/call/jsonrpc
    decorate with @services.jsonrpc the functions to expose
    supports xml, json, xmlrpc, jsonrpc, amfrpc, rss, csv
    """
    return service()


@auth.requires_signature()
def data():
    """
    http://..../[app]/default/data/tables
    http://..../[app]/default/data/create/[table]
    http://..../[app]/default/data/read/[table]/[id]
    http://..../[app]/default/data/update/[table]/[id]
    http://..../[app]/default/data/delete/[table]/[id]
    http://..../[app]/default/data/select/[table]
    http://..../[app]/default/data/search/[table]
    but URLs must be signed, i.e. linked with
      A('table',_href=URL('data/tables',user_signature=True))
    or with the signed load operator
      LOAD('default','data.load',args='tables',ajax=True,user_signature=True)
    """
    return dict(form=crud())

 

#currently saved or to previous version.
#Key bindings

#    Ctrl+S Save via Ajax
#    Ctrl+F11 Toggle Fullscreen
#    Ctrl-F / Cmd-F Start searching
#    Ctrl-G / Cmd-G Find Next
#    Shift-Ctrl-G / Shift-Cmd-G Find Previous
#    Shift-Ctrl-F / Cmd-Option-F Replace
#    Shift-Ctrl-R / Shift-Cmd-Option-F Replace All

#Powered by web2py™ created by Massimo Di Pierro ©2007-2013 - Admin language
