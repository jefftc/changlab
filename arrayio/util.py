"""

Functions:
num_headers  Guess the number of headers in a matrix.

"""

CHAR, INT, FLOAT, EMPTY = 1, 2, 4, 8
HEAD, SAMPLE, ANNOT, VALUE, BLANK = 1, 2, 4, 8, 16

def _rule_no_first_row_annots(matrix, num_rows, num_cols, datatype, semtype):
    # No ANNOT in the first row.
    changed = False
    for j in range(num_cols):
        if semtype[0][j] & ANNOT:
            semtype[0][j] ^= ANNOT
            changed = True
    return changed

def _rule_first_row_sample(matrix, num_rows, num_cols, datatype, semtype):
    # SAMPLE can only be in the first row.
    changed = False
    for i in range(1, num_rows):
        for j in range(num_cols):
            if semtype[i][j] & SAMPLE:
                semtype[i][j] ^= SAMPLE
                changed = True
    return changed

def _rule_first_row_col_head(matrix, num_rows, num_cols, datatype, semtype):
    # HEAD can only be in the first row or column.
    changed = False
    for i in range(1, num_rows):
        for j in range(1, num_cols):
            if semtype[i][j] & HEAD:
                semtype[i][j] ^= HEAD
                changed = True
    return changed

    ## # RULE 3: If the majority of the potential values in the matrix
    ## #         are floating points, then all values must be floating
    ## #         points.
    ## # This won't work.  E.g. if you use an algorithm to zero-fill
    ## # missing values.
    ## value_type_is = INT | FLOAT | EMPTY
    ## int_values = float_values = 0
    ## for i in range(num_rows):
    ##     for j in range(num_cols):
    ##         if not (semtype[i][j] & VALUE):
    ##             continue
    ##         if datatype[i][j] == INT:
    ##             int_values += 1
    ##         elif datatype[i][j] == FLOAT:
    ##             float_values += 1
    ## total = int_values + float_values
    ## if float_values >= total/2.0:
    ##     # Values must be FLOAT or EMPTY.
    ##     value_type_is = FLOAT | EMPTY
    ## for i in range(num_rows):
    ##     for j in range(num_cols):
    ##         if not (semtype[i][j] & VALUE):
    ##             continue
    ##         if value_type_is | datatype[i][j] != value_type_is:
    ##             semtype[i][j] ^= VALUE

def _rule_no_values_then_is_head(
    matrix, num_rows, num_cols, datatype, semtype):
    # If there are no VALUES in a column, then the first row, from
    # this down to the first column, must all be HEAD.
    changed = False
    for j in range(num_cols):
        any_values = False
        for i in range(num_rows):
            if semtype[i][j] & VALUE:
                any_values = True
                break
        if any_values:
            continue
        for jj in range(j+1):
            assert semtype[0][jj] & HEAD
            if semtype[0][jj] != HEAD:
                semtype[0][jj] = HEAD
                changed = True
    return changed

def _rule_no_broken_values(matrix, num_rows, num_cols, datatype, semtype):
    # In each row or column, the VALUEs can only appear at the end.
    changed = False
    for i in range(num_rows):
        in_value = True
        for j in range(num_cols-1, -1, -1):
            if in_value and not (semtype[i][j] & VALUE):
                in_value = False
            elif not in_value and (semtype[i][j] & VALUE):
                semtype[i][j] ^= VALUE
                changed = True
    for j in range(num_cols):
        in_value = True
        for i in range(num_rows-1, -1, -1):
            if in_value and not (semtype[i][j] & VALUE):
                in_value = False
            elif not in_value and (semtype[i][j] & VALUE):
                semtype[i][j] ^= VALUE
                changed = True
    return changed

def _rule_no_broken_head1(matrix, num_rows, num_cols, datatype, semtype):
    # In each row, the header must start from column 0.  There can't
    # be a cell with no HEAD followed by one with HEAD.  Same with
    # columns.
    changed = False
    for i in range(num_rows):
        in_header = True
        for j in range(num_cols):
            if in_header and not (semtype[i][j] & HEAD):
                in_header = False
            elif not in_header and (semtype[i][j] & HEAD):
                semtype[i][j] ^= HEAD
                changed = True
    for j in range(num_cols):
        in_header = True
        for i in range(num_rows):
            if in_header and not (semtype[i][j] & HEAD):
                in_header = False
            elif not in_header and (semtype[i][j] & HEAD):
                semtype[i][j] ^= HEAD
                changed = True
    return changed

def _rule_no_broken_head2(matrix, num_rows, num_cols, datatype, semtype):
    # If a cell is a HEAD, then all cells preceeding can only be HEAD.
    changed = False

    in_header = False
    for j in range(num_cols-1, -1, -1):
        if semtype[0][j] == HEAD:
            in_header = True
        elif in_header and semtype[0][j] != HEAD:
            semtype[0][j] = HEAD
            changed = True

    in_header = False
    for i in range(num_rows-1, -1, -1):
        if semtype[i][0] == HEAD:
            in_header = True
        elif in_header and semtype[i][0] != HEAD:
            semtype[i][0] = HEAD
            changed = True
    return changed


    ## # RULE 5: Label BLANK cells if there is a potential HEAD in the
    ## #         above it, a potential HEAD to its left, and no FLOATs
    ## #         anywhere to the right of it or below it.
    ## # This doesn't work:
    ## # <HEAD>  <HEAD>  <HEAD>  <SAMPLE>
    ## # <HEAD>  <BLANK> <BLANK> <VALUE>
    ## # <ANNOT> <ANNOT> <ANNOT> <VALUE>
    ## # If the <ANNOT> are numbers (e.g. GWEIGHT), then won't detect.
    ## could_be_BLANK = {}  # (i, j) -> 1
    ## for i in range(1, num_rows):
    ##     for j in range(1, num_rows):
    ##         if not (semtype[i][0] & HEAD) or not (semtype[0][j] & HEAD):
    ##             continue
    ##         any_floats = False
    ##         for ii in range(i+1, num_rows):
    ##             if datatype[ii][j] == FLOAT:
    ##                 any_floats = True
    ##                 break
    ##         if any_floats:
    ##             continue
    ##         for jj in range(j+1, num_cols):
    ##             if datatype[i][jj] == FLOAT:
    ##                 any_floats = True
    ##                 break
    ##         if any_floats:
    ##             continue
    ##         could_be_BLANK[(i, j)] = 1
    ## # Start with (1, 1) as BLANK.  Then add one row and column at a
    ## # time, making sure everything I added is blank.
    ## # X X X X
    ## # X X
    ## # X     X
    ## max_row = max_col = 1
    ## just_added_row = False
    ## while True:
    ##     if not just_added_row:
    ##         new_row, new_col = max_row+1, max_col
    ##         just_added_row = True
    ##     else:
    ##         new_row, new_col = max_row, max_col+1
    ##         just_added_row = False
        
    ##     all_blank = True
    ##     for i in range(1, new_row+1):
    ##         for j in range(1, new_col+1):
    ##             if (i, j) not in could_be_BLANK:
    ##                 all_blank = False
    ##     # If everything is BLANK, then accept the new rows and columns
    ##     # and try the next one.
    ##     if all_blank:
    ##         max_row, max_col = new_row, new_col
    ##         just_added_row = False
    ##     # If not everything is blank, and we just added a column, then
    ##     # we've already tried everything, and there's no more blanks.
    ##     elif not just_added_row:
    ##         break
    ## if (max_row, max_col) not in could_be_BLANK:
    ##     max_row = max_col = 0
    ## for i in range(1, max_row+1):
    ##     for j in range(1, max_col+1):
    ##         semtype[i][j] = BLANK
    ## for i in range(1, max_row+1):
    ##     semtype[i][0] = HEAD
    ## for j in range(1, max_col+1):
    ##     semtype[0][j] = HEAD
    ## for i in range(1, max_row+1):
    ##     for j in range(max_col+1, num_cols):
    ##         semtype[i][j] = ANNOT
    ## for j in range(1, max_col+1):
    ##     for i in range(max_row+1, num_rows):
    ##         semtype[i][j] = ANNOT

def _rule_no_broken_blank(matrix, num_rows, num_cols, datatype, semtype):
    # BLANKs can only be preceeded by BLANKs from (1, 1).  BLANKs must
    # have headers in the first row and column.
    changed = False
    blank_indexes = []  # list of (row, col)
    for i in range(num_rows):
        for j in range(num_cols):
            if semtype[i][j] & BLANK:
                blank_indexes.append((i, j))
    for i, j in blank_indexes:
        all_blank = True
        for ii in range(1, i):
            for jj in range(1, j):
                if not semtype[ii][jj] & BLANK:
                    all_blank = False
        if not all_blank:
            semtype[i][j] ^= BLANK
            changed = True
            continue
        if not semtype[i][0] & HEAD or not semtype[0][j] & HEAD:
            semtype[i][j] ^= BLANK
            changed = True
            continue
    return changed

def _rule_known_headers(matrix, num_rows, num_cols, datatype, semtype):
    # If the first row or column (except for (0, 0), because PCL files
    # allow different names) match known headers, then set them to
    # HEAD.

    KNOWN_COL_HEADERS = [
        "GID", "NA", "ID", "NAME", "LOCUSLINK",
        "GWEIGHT", "GORDER", "GCLUSTER"]
    KNOWN_ROW_HEADERS = ["GID", "AID", "EWEIGHT", "EORDER", "ACLUSTER"]

    changed = False

    if not num_rows or not num_cols:
        return changed
    if not semtype[0][0] & HEAD:
        return changed

    col_headers = [x.upper() for x in matrix[0]]
    for j in range(1, num_cols):
        if not semtype[0][j] & HEAD:
            break
        if col_headers[j] in KNOWN_COL_HEADERS:
            if semtype[0][j] != HEAD:
                semtype[0][j] = HEAD
                changed = True

    row_headers = [x[0].upper() for x in matrix]
    for i in range(1, num_rows):
        if not semtype[i][0] & HEAD:
            break
        if row_headers[i] in KNOWN_ROW_HEADERS:
            if semtype[i][0] != HEAD:
                semtype[i][0] = HEAD
                changed = True
    return changed

def _rule_no_values_by_head(
    matrix, num_rows, num_cols, datatype, semtype):
    # There are no VALUEs under column HEAD or to the right of row
    # HEAD.
    changed = False
    for j in range(num_cols):
        if semtype[0][j] != HEAD:
            break
        for i in range(num_rows):
            if semtype[i][j] & VALUE:
                semtype[i][j] ^= VALUE
                changed = True
    for i in range(num_rows):
        if semtype[i][0] != HEAD:
            break
        for j in range(num_cols):
            if semtype[i][j] & VALUE:
                semtype[i][j] ^= VALUE
                changed = True
    return changed

def _rule_head_around_blank(matrix, num_rows, num_cols, datatype, semtype):
    # RULE: If a cell has a HEAD on top and left, it must be BLANK.
    changed = False
    for i in range(1, num_rows):
        for j in range(1, num_cols):
            if not semtype[i][j] & BLANK:
                continue
            if semtype[i][j] == BLANK:
                continue
            if semtype[i][0] == HEAD and semtype[0][j] == HEAD:
                semtype[i][j] = BLANK
                changed = True
    return changed
    
def _rule_no_head_around_no_blank(
    matrix, num_rows, num_cols, datatype, semtype):
    # RULE: If a cell is not blank, then it cannot have a HEAD on the
    #       top and left.
    changed = False
    for i in range(1, num_rows):
        for j in range(1, num_cols):
            if semtype[i][j] & BLANK:
                continue
            assert semtype[i][0] != HEAD or semtype[0][j] != HEAD, \
                   "Ambiguous annotation."
            if semtype[i][0] == HEAD and semtype[0][j] & HEAD:
                semtype[0][j] ^= HEAD
                changed = True
            elif semtype[0][j] == HEAD and semtype[i][0] & HEAD:
                semtype[i][0] ^= HEAD
                changed = True
    return changed

def num_headers(matrix):
    """Return (# row headers, # col headers)."""
    # Try to find the number of rows and columns that contain header
    # information.

    # CASE 1:   No headers.  All <VALUES>
    # CASE 2:   1 row header, 1 column header.
    #           <HEAD>   <SAMPLE1>  <SAMPLE2>  [...]
    #           <ANNOT>   <VALUE>    <VALUE>
    # CASE 3:   1 row header, n column headers.
    #           <HEAD1>  <HEAD2>  <HEAD3>  <SAMPLE1>  <SAMPLE2>  <SAMPLE3>
    #           <ANNOT>  <ANNOT>  <ANNOT>   <VALUE>    <VALUE>    <VALUE>
    # CASE 4:   n row headers, 1 column headers.
    #           <HEAD1>  <SAMPLE1>  <SAMPLE2>  <SAMPLE3>
    #           <HEAD4>   <ANNOT>    <ANNOT>    <ANNOT>
    #           <HEAD5>   <ANNOT>    <ANNOT>    <ANNOT>
    #           <ANNOT>   <VALUE>    <VALUE>    <VALUE>
    # CASE 5:   n row headers, n column headers.
    #           <HEAD1>  <HEAD2>  <HEAD3>  <SAMPLE1>  <SAMPLE2>  <SAMPLE3>
    #           <HEAD4>  <BLANK>  <BLANK>   <ANNOT>    <ANNOT>    <ANNOT>
    #           <HEAD5>  <BLANK>  <BLANK>   <ANNOT>    <ANNOT>    <ANNOT>
    #           <ANNOT>  <ANNOT>  <ANNOT>   <VALUE>    <VALUE>    <VALUE>
    #
    #               1    2     4      8
    #  1 HEAD     CHAR  INT  FLOAT
    #  2 SAMPLE   CHAR  INT  FLOAT
    #  4 ANNOT    CHAR  INT  FLOAT  EMPTY
    #  8 VALUE          INT  FLOAT  EMPTY
    # 16 BLANK                      EMPTY
    #
    # Challenges:
    # - It's hard to distinguish between ANNOT, VALUE, and BLANK when
    #   they are EMPTY.
    # - It's hard to distinguish between HEADs and SAMPLEs.
    #
    # RULE: No ANNOT in the first row.
    # RULE: SAMPLE can only be in the first row.
    # RULE: HEAD can only be in the first row or column.
    # RULE: If there are no VALUES in a column, then the first row,
    #       from this down to the first column, must all be HEAD.
    # RULE: In each row, the header must start from column 0.  There
    #       can't be a cell with no HEAD followed by one with HEAD.
    # RULE: If a cell is a HEAD, then all cells preceeding can only be
    #       HEAD.
    # RULE: BLANKs can only be preceeded by BLANKs from (1, 1).
    #       BLANKs must have headers in the first row and column.
    # RULE: In each row or column, the VALUEs can only appear from
    #       the end.
    # RULE: If the first row or column (except for (0, 0), because
    #       PCL files allow different names) match known headers,
    #       then set them to HEAD.
    # RULE: There are no VALUEs under column HEAD or to the right of
    #       row HEAD.
    # RULE: If a cell has a HEAD on top and left, that cell must be
    #       BLANK.
    # RULE: If a cell is not blank, then it cannot have a HEAD on the
    #       top and left.
    RULES = [
        _rule_no_first_row_annots,
        _rule_first_row_sample,
        _rule_first_row_col_head,
        _rule_no_values_then_is_head,
        _rule_no_broken_head1,
        _rule_no_broken_head2,
        _rule_no_broken_blank,
        _rule_no_broken_values,
        _rule_known_headers,
        _rule_no_values_by_head,
        _rule_head_around_blank,
        _rule_no_head_around_no_blank,
        ]

    if not matrix:
        return 0, 0
    num_rows, num_cols = len(matrix), len(matrix[0])
    # Make sure each row contains the same number of columns.
    for row in matrix:
        assert len(row) == num_cols, "matrix row length mismatch"

    # Figure out the data type for each cell in the matrix.
    CHAR, INT, FLOAT, EMPTY = 1, 2, 4, 8
    datatype = [[None]*num_cols for i in range(num_rows)]
    for i in range(num_rows):
        for j in range(num_cols):
            x = matrix[i][j]
            if x.strip() == "":
                dt = EMPTY
            elif _is_int(x):
                dt = INT
            elif _is_float(x):
                dt = FLOAT
            else:
                dt = CHAR
            datatype[i][j] = dt

    # Make an initial guess at the semantic types of each cell.
    HEAD, SAMPLE, ANNOT, VALUE, BLANK = 1, 2, 4, 8, 16
    semtype = [[0]*num_cols for i in range(num_rows)]
    for i in range(num_rows):
        for j in range(num_cols):
            x = datatype[i][j]
            if x == CHAR:
                st = HEAD | SAMPLE | ANNOT
            elif x == INT:
                st = HEAD | SAMPLE | ANNOT | VALUE
            elif x == FLOAT:
                st = HEAD | SAMPLE | ANNOT | VALUE
            elif x == EMPTY:
                st = ANNOT | VALUE | BLANK
            else:
                raise AssertionError
            semtype[i][j] = st

    # Apply the rules to guess the right types of each cell of the
    # matrix.
    changed = True
    while changed:
        changed = False
        for rule_fn in RULES:
            c = rule_fn(matrix, num_rows, num_cols, datatype, semtype)
            changed = changed or c

    # Look for the VALUEs.  Start looking at the bottom right of the
    # MATRIX, and add one column and row at a time.
    first_row, first_col = num_rows-1, num_cols-1
    just_added_row = False
    while True:
        if not just_added_row:
            new_row, new_col = first_row-1, first_col
            just_added_row = True
        else:
            new_row, new_col = first_row, first_col-1
            just_added_row = False
        # Make sure the rows and cols are in bounds.
        if new_row < 0 or new_col < 0:
            if just_added_row:
                continue
            break

        all_values = True
        for i in range(new_row, num_rows):
            for j in range(new_col, num_cols):
                if not semtype[i][j] & VALUE:
                    all_values = False
        # If everything is a VALUE, then accept the new rows and
        # columns and try the next one.
        if all_values:
            first_row, first_col = new_row, new_col
            just_added_row = False
        # If not everything is a value, and we just added a column,
        # then we've already tried everything, and there's no more
        # values.
        elif not just_added_row:
            break
    if not semtype[first_row][first_col] & VALUE:
        # There are no values.
        first_row, first_col = num_rows, num_cols
    hrows, hcols = first_row, first_col

    #print datatype
    #print semtype
    return hrows, hcols

def _all_numeric(vec):
    for n in vec:
        if not _is_numeric(n):
            return False
    return True

def _is_numeric(n):
    # empty strings are not numeric.
    if n == "":
        return False
    try:
        float(n)
    except ValueError, x:
        return False
    return True

def _is_int(n):
    try:
        int(n)
    except ValueError, x:
        return False
    return True

def _is_float(n):
    try:
        float(n)
    except ValueError, x:
        return False
    return True

def _is_float_not_int(n):
    if _is_int(n):
        return False
    try:
        float(n)
    except ValueError, x:
        return False
    return True
