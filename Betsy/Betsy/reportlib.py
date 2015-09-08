# For help with generating HTML report files.

def extract_filenames(antecedents, out_path, rename_outfile_fn=None):
    # Return list of (in_file, out_file, in_filename, out_filename)
    # where the filenames are the full paths, while the files are just
    # the file part.
    import os
    
    filenames = []
    for i, data_node in enumerate(antecedents):
        in_filename = data_node.identifier
        in_path, in_file = os.path.split(in_filename)
        out_file = in_file
        if rename_outfile_fn:
            out_file = rename_outfile_fn(i, in_file)
        out_filename = os.path.join(out_path, out_file)
        x = in_file, out_file, in_filename, out_filename
        filenames.append(x)
    return filenames


def copy_file_or_path(in_filename, out_filename):
    import os
    import shutil

    if os.path.exists(out_filename):
        # Need to clean up the path so that old files aren't left
        # over.
        if os.path.isdir(out_filename):
            shutil.rmtree(out_filename)
        else:
            os.unlink(out_filename)
    
    if os.path.isdir(in_filename):
        shutil.copytree(in_filename, out_filename)
    else:
        shutil.copyfile(in_filename, out_filename)
