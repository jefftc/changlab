"""

Functions:
read_cls_file
write_cls_file

resolve_classes

"""

def read_cls_file(filename):
    # Return tuple (class_names, classes).  class_names is a list
    # containing the names of the classes.  classes is a list of the
    # classes given in the file.  Classes are either from class_names,
    # or an integer from [0, class_names).
    # 
    # Limitations:
    # Only handles categorical CLS files with 2 classes.
    import filelib
    
    # Space or tab-delimited format.
    # <num samples> <num classes> 1
    # # <class name 0> <class name 1> ...
    # <0/1 or class name> ...
    handle = filelib.openfh(filename)
    x = [x for x in handle if x.strip()]
    assert len(x) == 3, "CLS file should contain 3 lines."
    line1, line2, line3 = x

    # Parse the first line.
    x = line1.strip().split()
    assert len(x) == 3
    assert x[2] == "1"
    num_samples, num_classes = int(x[0]), int(x[1])

    # Parse the second line.
    x = line2.strip().split()
    assert x
    assert x[0] == "#"
    assert len(x) == num_classes+1
    class_names = x[1:]

    # Parse the third line.
    x = line3.strip().split()
    assert len(x) == num_samples
    classes = x
    for i, x in enumerate(classes):
        if x in class_names:
            continue
        try:
            x = int(x)
        except ValueError:
            assert False, "Invalid class: %s" % x
        assert x >= 0 and x < num_classes
        classes[i] = x
    return class_names, classes
    
def write_cls_file(outhandle, name0, name1, classes):
    # Only handles categorical CLS files with 2 classes.
    # classes should be a list of 0/1 or class names.
    from genomicode import hashlib
    
    # Check the classes variable.
    assert classes
    #for x in classes:
    #    assert x in [0, 1, "0", "1", name0, name1]
    uniq_classes = []
    for x in classes:
        if x not in uniq_classes:
            uniq_classes.append(x)
    assert len(uniq_classes) == 2, "Need exactly 2 classes."
    sorted_classes = sorted(map(str, uniq_classes))
    assert sorted_classes in [["0", "1"], sorted([name0, name1])]
    # Make sure order of the classes is consistent with the names.
    assert str(uniq_classes[0]) in ["0", name0], "classes out of order"
    
    if type(outhandle) is type(""):
        outhandle = open(outhandle, 'w')

    # Space or tab-delimited format.
    # <num samples> <num classes> 1
    # # <class name 0> <class name 1> ...
    # <0/1 or class name> ...
    num_samples = len(classes)
    x = [num_samples, 2, 1] + [""]*(num_samples-3)
    print >>outhandle, "\t".join(map(str, x))

    hname0, hname1 = hashlib.hash_var(name0), hashlib.hash_var(name1)
    assert hname0 != hname1
    x = ["#", hname0, hname1] + [""]*(num_samples-3)
    print >>outhandle, "\t".join(map(str, x))

    print >>outhandle, "\t".join(map(str, classes))


def resolve_classes(MATRIX, indexes1, indexes2, count_headers, name1, name2):
    # indexes1 is a string.  indexes2 is a string or None.
    # Return name1, name2, classes.  classes is 0, 1, or None.
    from genomicode import parselib
    
    max_index = MATRIX.ncol()
    num_headers = len(MATRIX._row_names)
    assert max_index, "empty matrix"
    
    assert indexes1 and type(indexes1) is type("")

    I1 = []
    for s, e in parselib.parse_ranges(indexes1):
        if count_headers:
            s, e = s - num_headers, e - num_headers
        assert s >= 1, "Index out of range: %s" % s
        assert e <= max_index, "Index out of range: %s" % e
        s, e = s - 1, min(e, max_index)
        I1.extend(range(s, e))
    
    I2 = []
    if indexes2:
        for s, e in parselib.parse_ranges(indexes2):
            if count_headers:
                s, e = s - num_headers, e - num_headers
            assert s >= 1, "Index out of range: %s" % s
            assert e <= max_index, "Index out of range: %s" % e
            s, e = s - 1, min(e, max_index)
            I2.extend(range(s, e))
    else:
        # If indexes2 not given, then I2 should be every index that's
        # not in I1.
        I2 = [i for i in range(max_index) if i not in I1]
        
    # Make sure no overlap between I1 and I2.
    for i in I1:
        assert i not in I2, "Overlap in classes."

    # Provide default group names.
    # If there is only 1 index, then use the sample name from the
    # matrix.
    col_header = None
    if MATRIX.col_names():
        col_header = MATRIX.col_names()[0]
    if not name1 and col_header and len(I1) == 1:
        name1 = MATRIX.col_names(col_header)[I1[0]]
    if not name2 and col_header and len(I2) == 1:
        name2 = MATRIX.col_names(col_header)[I2[0]]
    
    name1 = name1 or "group1"
    name2 = name2 or "group2"
    if name1 == name2:
        name1 = "%s-1" % name1
        name2 = "%s-2" % name2

    classes = [None]*MATRIX.ncol()
    for i in I1:
        classes[i] = 0
    for i in I2:
        classes[i] = 1

    x = name1, name2, classes
    return x
