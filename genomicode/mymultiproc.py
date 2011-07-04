"""

Functions:
run

"""



def run(jobs, num_procs=4, lock_keyword=None):
    # Each job is a tuple of (function, args, keywds).
    # If lock_keyword is provided, then I will make a lock and pass it
    # to the function with this name.
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
                # Should raise exception here instead?
                print >>sys.stderr, "ERROR: %s" % str(x)
                continue
            results[i] = x
        if all_done:
            break
        time.sleep(0.01)

    return results
