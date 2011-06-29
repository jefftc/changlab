"""

Functions:
run

"""



def run(jobs, num_procs=4, lock_keyword=None):
    # Each job is a tuple of (function, args, keywds).
    # If lock_keyword is provided, then I will make a lock and pass it
    # to the function with this name.
    import multiprocessing

    assert num_procs >= 1 and num_procs <= 1024
    
    manager = multiprocessing.Manager()
    pool = multiprocessing.Pool(num_procs)

    lock = None
    if lock_keyword is not None:
        lock = manager.Lock()
        
    results = []
    for i, job in enumerate(jobs):
        assert len(job) == 3
        fn, args, keywds = job

        if lock_keyword is not None:
            keywds = keywds.copy()
            keywds[lock_keyword] = lock
        
        if num_procs == 1:
            fn(*args, **keywds)
        else:
            x = pool.apply_async(fn, args, keywds)
            results.append(x)
    pool.close()
    pool.join()
    for x in results:
        x.get()
