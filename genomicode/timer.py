"""

Functions:
wait

"""

TIMERS = {}   # name -> time

def wait(delay, name=None):
    global TIMERS
    import time

    if delay is None:
        delay = 2
    if name is None:
        name = "default"

    how_long = TIMERS.get(name, 0) + delay - time.time()
    if how_long > 0:
        time.sleep(how_long)
    TIMERS[name] = time.time()
