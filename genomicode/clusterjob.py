import os, sys

class ClusterJob:
    def load(self):
        # Should estimate the computational load of this job in
        # arbitrary units.
        raise NotImplementedError
    def join(self, job):
        # Should return the result when this job is joined with another.
        raise NotImplementedError
    def split(self):
        # Should split this job in half and return a tuple of job1,
        # job2.
        raise NotImplementedError


def _find_smallest_job(jobs):
    # Return the index of the smallest job.
    min_load = min_i = None
    for i, job in enumerate(jobs):
        load = job.load()
        if min_load is None or load < min_load:
            min_load, min_i = load, i
    return min_i

def _find_biggest_job(jobs):
    # Return the index of the biggest job.
    max_load = max_i = None
    for i, job in enumerate(jobs):
        load = job.load()
        if max_load is None or load > max_load:
            max_load, max_i = load, i
    return max_i
    
def balance(all_jobs, num_jobs=None):
    # Return a list of jobs, roughly balanced.
    assert len(all_jobs), "no jobs specified"
    if not num_jobs:
        num_jobs = len(all_jobs)

    # Balanced the jobs using simple heuristics.
    jobs = all_jobs[:]
    zzz = 0
    zzz_stop = num_jobs + max(num_jobs, len(jobs))
    while zzz < zzz_stop or len(jobs) != num_jobs:
        zzz += 1
        print "%d:%d %d %d" % (zzz, zzz_stop, len(jobs), num_jobs)
        sys.stdout.flush()
        if len(jobs) < num_jobs:
            i = _find_biggest_job(jobs)
            job1, job2 = jobs[i].split()
            jobs[i] = job1
            jobs.append(job2)
            print "Split %d" % i
        elif len(jobs) > num_jobs:
            import random
            #i = _find_smallest_job(jobs)
            i = random.randint(0, len(jobs)-1)
            job = jobs.pop(i)
            #j = _find_smallest_job(jobs)
            j = random.randint(0, len(jobs)-1)
            jobs[j] = jobs[j].join(job)
            print "Merge %d %d" % (i, j)
        else:
            i = _find_smallest_job(jobs)
            j = _find_biggest_job(jobs)
            if i == j:
                continue
            job = jobs[i].join(jobs[j])
            job1, job2 = job.split()
            jobs[i], jobs[j] = job1, job2
            print "Balance %d %d" % (i, j)
    return jobs

def format_dscr(name, path, high_priority, command):
    priority = "#$ -pe low-all 1"
    if high_priority:
        priority = "#$ -pe high 1\n#$ -l highprio"

    lines = []
    w = lines.append
    w("#!/bin/tcsh\n")
    w("#\n")
    w("#$ -S /bin/tcsh -cwd\n")
    w("#$ -o %s/%s.out -j y\n" % (path, name))
    w("%s\n" % priority)
    w("cd %s\n" % path)
    w("%s\n" % command)
    return "".join(lines)
