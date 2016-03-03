"""Run commands in parallel.

Methods:
sshell  Run a single shell command.
pshell  Run a bunch of shell commands in parallel.
pyfun   Run python functions in parallel.

quote

"""
    
def sshell(command, path=None, ignore_nonzero_exit=False):
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
        x = _sshell_h(command, ignore_nonzero_exit)
    finally:
        os.chdir(cwd)
    return x
    

def _sshell_h(command, ignore_nonzero_exit):
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
        if ignore_nonzero_exit:
            x = x3
        else:
            msg = "\n".join([x1, x2, x3])
            raise AssertionError, msg
    return x


def pshell(commands, max_procs=None, path=None):
    # commands is a list of shell commands to run.  Return the output
    # as a single string.
    import os

    assert type(commands) is not type("")

    cwd = os.getcwd()
    try:
        if path is not None:
            if not os.path.exists(path):
                os.mkdir(path)
            os.chdir(path)
        x = _pshell_h(commands, max_procs)
    finally:
        os.chdir(cwd)
    return x
    

def _pshell_h(commands, max_procs):
    import subprocess
    
    #cat run01.sh | parallel -j <num_cores>
    parallel_bin = _find_parallel()
    assert parallel_bin, "Not found: parallel"
    cmd = [
        parallel_bin,
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


def pyfun(jobs, num_procs=4, lock_keyword=None, DELAY=0.1,
        propogate_exception=True):
    # Each job is a tuple of (function, args, keywds).
    # args and keywds can be None.
    # If lock_keyword is provided, then I will make a lock and pass it
    # to the function with this name.
    # DELAY is the number of seconds to wait before polling for
    # completed processes.
    import sys
    import time
    import multiprocessing

    assert num_procs >= 1 and num_procs <= 1024
    
    manager = multiprocessing.Manager()
    pool = multiprocessing.Pool(num_procs)

    lock = None
    if lock_keyword is not None:
        lock = manager.Lock()

    results = [None] * len(jobs)
    procs = []
    for i, job in enumerate(jobs):
        assert len(job) == 3
        fn, args, keywds = job

        if args is None:
            args = ()
        if keywds is None:
            keywds = {}

        if lock_keyword is not None:
            keywds = keywds.copy()
            keywds[lock_keyword] = lock
        
        if num_procs == 1:
            x = fn(*args, **keywds)
            results[i] = x
        else:
            x = pool.apply_async(fn, args, keywds)
            procs.append(x)
    pool.close()
    # Need this, otherwise get exception:
    # Exception RuntimeError: RuntimeError('cannot join current
    # thread',) in <Finalize object, dead> ignored
    pool.join()
    assert len(procs) == 0 or len(procs) == len(jobs)

    done = [False] * len(procs)
    while 1:
        all_done = True
        for i in range(len(done)):
            if done[i]:
                continue
            if not procs[i].ready():
                all_done = False
                continue
            done[i] = True
            try:
                x = procs[i].get()
            except (SystemError, KeyboardInterrupt, MemoryError), x:
                raise
            except Exception, x:
                if propogate_exception:
                    raise
                # Should raise exception here instead?
                print >>sys.stderr, "ERROR: %s" % str(x)
                continue
            results[i] = x
        if all_done:
            break
        time.sleep(DELAY)

    return results


def quote(s, always_quote=False):
    BAD_CHARS = " \\"

    if type(s) in [type(0), type(0.0)]:
        s = str(s)
    needs_quote = False
    if not always_quote:
        for x in BAD_CHARS:
            if x in s:
                needs_quote = True
                break
    if always_quote or needs_quote:
        s = "'" + s.replace("'", "'\\''") + "'"
    return s



def _find_parallel():
    # Return the full path for GNU parallel.
    import os
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
