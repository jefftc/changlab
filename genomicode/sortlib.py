def hash_natural(text):
    import re
    def convert(x):
        if x.isdigit():
            return int(x)
        return x
    # This is not compatible with Python 2.3.4.
    #convert = lambda x: int(x) if x.isdigit() else x
    return [convert(x) for x in re.split('([0-9]+)', text)]

def sort_natural(L):
    schwartz = [(hash_natural(x), x) for x in L]
    schwartz.sort()
    x = [x[-1] for x in schwartz]
    return x
