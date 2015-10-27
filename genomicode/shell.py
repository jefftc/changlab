"""Run shell commands.

Methods:
single    Run a single shell command.
parallel  Run a bunch of shell commands in parallel.

quote

"""

def single(command, path=None):
    # command is a string or list (see subprocess) for a single
    # command.  If path is not None, will run in this path and return
    # to the current working directory when done.  Return the output
    # as a single string.
    import os

    cwd = os.getcwd()
    try:
        if path is not None:
            if not os.path.exists(path):
                os.mkdir(path)
            os.chdir(path)
        x = _single_h(command)
    finally:
        os.chdir(cwd)
    return x
    

def _single_h(command):
    import subprocess
    if type(command) != type(""):
        command = " ".join(command)
    try:
        x = subprocess.check_output(
            command, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError, x:
        x1 = "Non-zero exit status [%d]:" % x.returncode
        x2 = command
        x3 = x.output
        msg = "\n".join([x1, x2, x3])
        raise AssertionError, msg
    return x


def _parallel_h(commands, max_procs):
    import subprocess
    
    #cat run01.sh | parallel -j <num_cores>
    cmd = [
        "parallel",
        ]
    if max_procs is not None:
        assert max_procs >= 0 and max_procs < 100
        cmd.extend(["-j", str(max_procs)])

    p = subprocess.Popen(
        cmd, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    x = "\n".join(commands)
    x = p.communicate(x)
    data_out, data_err = x

    x = "%s\n%s" % (data_out, data_err)
    x = x.strip()
    return x


def parallel(commands, max_procs=None, path=None):
    # commands is a list of shell commands to run.  Return the output
    # as a single string.
    import os

    cwd = os.getcwd()
    try:
        if path is not None:
            if not os.path.exists(path):
                os.mkdir(path)
            os.chdir(path)
        x = _parallel_h(commands, max_procs)
    finally:
        os.chdir(cwd)
    return x
    


def quote(s):
    return "'" + s.replace("'", "'\\''") + "'"


