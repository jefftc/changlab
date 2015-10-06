"""

Functions:
write_scatterplot

"""

def write_scatterplot(filename, DATA, rownames=None):
    # DATA is a dictionary where the keys are the NAMEs of each series
    # and the values are a two-dimensional list of xy coordinates.
    # rownames is a dictionary where the keys are the NAMEs and the
    # values are lists of names (same length as the matrices).

    # Prism matrix format has one column of x's, and one column of Ys
    # for each series.
    # X  Y1  Y2  ...

    matrix = []
    series_all = sorted(DATA)
    for i, name in enumerate(series_all):
        xy = DATA[name]
        mat = [[""]*(len(series_all)+1) for j in range(len(xy))]
        for j in range(len(xy)):
            x, y = xy[j]
            mat[j][0] = x
            mat[j][i+1] = y
        matrix.extend(mat)


    # If there are row names, then add it to the matrix.
    if rownames:
        all_names = []
        for name in series_all:
            all_names.extend(rownames[name])
        assert len(all_names) == len(matrix)
        matrix_new = []
        for i in range(len(matrix)):
            x = [all_names[i]] + matrix[i]
            matrix_new.append(x)
        matrix = matrix_new

    # Add the header.
    header = ["X"] + series_all
    if rownames:
        header = ["Title"] + header
    matrix = [header] + matrix

    # Write out the matrix.
    handle = open(filename, 'w')
    for x in matrix:
        print >>handle, "\t".join(map(str, x))

