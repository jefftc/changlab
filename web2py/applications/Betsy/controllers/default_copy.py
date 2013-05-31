# coding: utf8
## -*- coding: utf-8 -*-

upload_path = '/home/xchen/chencode/upload/'
def status():
    return dict(request=request,session=session,response=response)

   
@auth.requires_login()
def index():
    form1 = FORM(INPUT(_type='submit',_value='normalize file',_name='normalize'))
    form2 = FORM(INPUT(_type='submit',_value='cluster genes',_name='cluster'))
    form3 = FORM(INPUT(_type='submit',_value='make heatmap',_name='heatmap'))
    form4 = FORM(INPUT(_type='submit',_value='classification',_name='classification'))
    form5 = FORM(INPUT(_type='submit',_value='batch_effect_remove',_name='batch'))
    form6 = FORM(INPUT(_type='submit',_value='differential expressed gene analysis',_name='diffgene'))
    form7 = FORM(INPUT(_type='submit',_value='geneset analysis',_name='geneset'))
    form_queue = FORM(INPUT(_type='submit',_value='show jobs status',_name='status'))
    if form_queue.process(formname='status').accepted:
        redirect(URL('show_status'))
    if form1.process(formname='normalize').accepted:
        redirect(URL('normalize_file'))
    if form2.process(formname='cluster').accepted:
        redirect(URL('cluster_genes'))
    if form3.process(formname='heatmap').accepted:
        redirect(URL('heatmap'))
    if form4.process(formname='classification').accepted:
        redirect(URL('classification'))
    if form5.process(formname='batch').accepted:
        redirect(URL('batch_effect'))
    if form6.process(formname='diffgene').accepted:
        redirect(URL('diffgene'))
    if form7.process(formname='geneset').accepted:
        redirect(URL('geneset'))
    return dict(form1=form1,form2=form2,form3=form3,form4=form4,form5=form5,form6=form6,form7=form7,form_queue=form_queue)


def normalize_file():
    fields = create_form(normalize_inputs,normalize_parameters,normalize_default)
    form=SQLFORM.factory(*fields)
    form_queue = FORM(INPUT(_type='submit',_value='show jobs status'))
    if form.process().accepted:
        form_value = form.vars
        input_file,command_parameters = create_command(form_value)
        session.parameters = command_parameters
        session.input_file = input_file
        session.protocol = 'normalize_file'
        redirect(URL('submit_job'))
    if form_queue.process().accepted:
        redirect(URL('show_status'))
    return dict(form=form,form_queue=form_queue)
     

def cluster_genes():
    fields = create_form(cluster_inputs,cluster_parameters,cluster_default)
    form=SQLFORM.factory(*fields)
    form_queue = FORM(INPUT(_type='submit',_value='show jobs status'))
    if form.process().accepted:
        form_value = form.vars
        input_file,command_parameters = create_command(form_value)
        session.parameters = command_parameters
        session.input_file = input_file
        session.protocol = 'cluster_genes'
        redirect(URL('submit_job'))
    if form_queue.process().accepted:
        redirect(URL('show_status'))
    return dict(form=form,form_queue=form_queue)

 
def heatmap():
    fields = create_form(heatmap_inputs,heatmap_parameters,heatmap_default)
    form=SQLFORM.factory(*fields)
    form_queue = FORM(INPUT(_type='submit',_value='show jobs status'))
    if form.process().accepted:
        form_value = form.vars
        input_file,command_parameters = create_command(form_value)
        session.parameters = command_parameters
        session.input_file = input_file
        session.protocol = 'make_heatmap'
        redirect(URL('submit_job'))
    if form_queue.process().accepted:
        redirect(URL('show_status'))
    return dict(form=form,form_queue=form_queue)
   

def classification():
    fields = create_form(classification_inputs,classification_parameters,classification_default)
    form=SQLFORM.factory(*fields)
    form_queue = FORM(INPUT(_type='submit',_value='show jobs status'))
    if form.process().accepted:
        form_value = form.vars
        input_file,command_parameters = create_command(form_value)
        session.parameters = command_parameters
        session.input_file = input_file
        session.protocol = 'classification'
        redirect(URL('submit_job'))
    if form_queue.process().accepted:
        redirect(URL('show_status'))
    return dict(form=form,form_queue=form_queue)

 
def batch_effect():
    fields = create_form(batch_inputs,batch_parameters,batch_default)
    form=SQLFORM.factory(*fields)
    form_queue = FORM(INPUT(_type='submit',_value='show jobs status'))
    if form.process().accepted:
        form_value = form.vars
        input_file,command_parameters = create_command(form_value)
        session.parameters = command_parameters
        session.input_file = input_file
        session.protocol = 'batch_effect_remove'
        redirect(URL('submit_job'))
    if form_queue.process().accepted:
        redirect(URL('show_status'))
    return dict(form=form,form_queue=form_queue)

    
def diffgene():
    fields = create_form(diffgene_inputs,diffgene_parameters,diffgene_default)
    form=SQLFORM.factory(*fields)
    form_queue = FORM(INPUT(_type='submit',_value='show jobs status'))
    if form.process().accepted:
        form_value = form.vars
        input_file,command_parameters = create_command(form_value)
        session.parameters = command_parameters
        session.input_file = input_file
        session.protocol = 'differential_expressed_gene_analysis'
        redirect(URL('submit_job'))
    if form_queue.process().accepted:
        redirect(URL('show_status'))
    return dict(form=form,form_queue=form_queue)

    
def geneset():
    fields = create_form(geneset_inputs,geneset_parameters,geneset_default)
    form=SQLFORM.factory(*fields)
    form_queue = FORM(INPUT(_type='submit',_value='show jobs status'))
    if form.process().accepted:
        form_value = form.vars
        input_file,command_parameters = create_command(form_value)
        session.parameters = command_parameters
        session.input_file = input_file
        session.protocol = 'geneset_analysis'
        redirect(URL('submit_job'))
    if form_queue.process().accepted:
        redirect(URL('show_status'))
    return dict(form=form,form_queue=form_queue)

    

def show_status():
    import subprocess
    all_cmd = ['python','/home/xchen/queue/queue.py','--list']    
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
                Ended=lines[4],Submitted=lines[5],Command=lines[6])
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
        pro_file = os.path.join('/home/xchen/queue/processing_info/', hash_string+'.txt')
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
    cmd = ['python', '/home/xchen/chencode/Betsy/scripts/run_protocol.py', 
           '--protocol', session.protocol] + session.input_file  + session.parameters
    all_cmd = ['python','/home/xchen/queue/queue.py' ,' ' .join(cmd)]
    process = subprocess.Popen(all_cmd, shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)    
    process.wait()                      
    redirect(URL('show_status'))
    return dict()

    
def result_folder():
    path = os.path.join('/home/xchen/chencode/examples/sample_analysis',request.args[0]+'/'+request.args[1])
    return open(path,'rb').read()

    
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
