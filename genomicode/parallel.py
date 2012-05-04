# Functions for running GNU parallel.

# Should deprecate this in favor of multiprocessing module.

import os, sys

def find():
    # Return the full path for GNU parallel.
    import subprocess
    
    PATHS = ["/usr/local/bin", "/usr/bin", "/bin",
             "/opt/local/bin", "/opt/bin"]
    parallel_bin = "parallel"

    for path in PATHS:
        filename = os.path.join(path, parallel_bin)

        # Check if the path exists.
        if not os.path.exists(filename):
            continue
        # Check if this is GNU parallel.
        cmd = "%s --version" % filename
        p = subprocess.Popen(
            cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        w, r = p.stdin, p.stdout
        w.close()
        output = r.read()
        if not output.startswith("GNU parallel"):
            continue
        return filename

    return None

# Problem: How to detect if the process failed?
def run_commands(commands, num_procs):
    # commands is a list of shell commands (as strings).  num_procs is
    # the number of processes to run.
    import subprocess
    
    assert num_procs > 0 and num_procs < 128

    parallel_bin = find()
    assert parallel_bin, "parallel not found"
    
    cmd = "%s -j %d" % (parallel_bin, num_procs)
    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    w, r = p.stdin, p.stdout

    for cmd in commands:
        w.write("%s\n" % cmd)
    w.close()

    output = r.read()
    return output
