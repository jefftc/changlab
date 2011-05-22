# Functions that deal with compressed archives.
#
# Functions:
# zip_path
# 
# unzip_list   Return list of (shortname, fullname).
# unzip_dict   Return dict of shortname -> fullname.

import os, sys

def zip_path(file_or_path, zipfile=None, noclobber=True):
    import stat
    import subprocess
    
    assert os.path.exists(file_or_path), "Missing path: %s" % file_or_path
    root, name  = os.path.split(file_or_path)

    # Figure out a good name for the zipfile, if not given.
    if not zipfile:
        zipfile = "%s.zip" % os.path.splitext(name)[0]
    assert not (noclobber and os.path.exists(zipfile)), "Zipfile exists."
    
    # Find "zip".  Somehow the "zip" on the path for cagt does not
    # find the right command.  Will get errors:
    # /bin/sh: zip: command not found
    # Check for "zip" in the usual places.  If it doesn't exist,
    # then try without a path.
    ZIP_PATHS = ["/usr/local/bin/zip", "/usr/bin/zip", "/bin/zip"]
    zipcmd = "zip"
    for path in ZIP_PATHS:
        if os.path.exists(path):
            zipcmd = path
            break

    cwd = os.getcwd()
    try:
        if root:
            os.chdir(root)

        cmd = "%s -r '%s' '%s'" % (zipcmd, zipfile, name)
        #w, r = os.popen4(cmd)
        #print cmd
        p = subprocess.Popen(
            cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        w, r = p.stdin, p.stdout
        w.close()

        # Check for errors in the zip output.
        output = r.read()
        assert output.find("command not found") < 0, x

        # Make sure the zip file exists.
        # This doesn't work.  exists_nz will also check for the file
        # without ".zip" and will find it.
        #assert filelib.exists_nz(zipfile), "Failed to archive trash."
        assert os.path.exists(zipfile) and os.stat(zipfile)[stat.ST_SIZE] > 0,\
               "Failed to make archive."

        # Clear out the old directory.
        if zipfile != name:
            os.system("rm -rf '%s'" % name)
    finally:
        os.chdir(cwd)

    return output

def unzip_list(filename):
    import zipfile

    assert os.path.exists(filename), "I could not find file %s." % filename
    assert zipfile.is_zipfile(filename)

    zfile = zipfile.ZipFile(filename)
    fullnames = zfile.namelist()
    shortnames = [os.path.split(x)[1] for x in fullnames]
    return zip(shortnames, fullnames)

def unzip_dict(filename):
    s2f = {}
    for x in unzip_list(filename):
        shortname, fullname = x
        assert shortname not in s2f, "Duplicate short name: %s" % shortnam
        s2f[shortname] = fullname
    return s2f
