"""

Functions:
run

"""

import os

def run(script, matlab_bin=None, working_path=None):
    import subprocess
    import config
    
    matlab_bin = matlab_bin or config.matlab or "matlab"
    matlab_args = ["-nosplash", "-nodesktop", "-nodisplay"]
    x = " ".join(matlab_args)
    matlab_cmd = "%s %s" % (matlab_bin, x)

    cwd = os.getcwd()
    try:
        if working_path:
            os.chdir(working_path)
        p = subprocess.Popen(
            matlab_cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        w, r = p.stdin, p.stdout
        w.writelines(script)
        w.close()
        output = r.read()
    finally:
        if working_path:
            os.chdir(cwd)
    return output
