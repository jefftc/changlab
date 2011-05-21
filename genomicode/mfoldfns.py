"""

Functions:
calc_free_energy
mfold

"""
import os, sys
from filefns import openfh

def calc_free_energy(sequence, NA=None, T=None):
    # Returns the delta G of the best folding.  If no foldings, then
    # returns 0.
    import time
    import tempfile
    import genomefns

    tempdir = tempfile.mkdtemp(dir=".")
    assert os.path.exists(tempdir)

    # mfold generates temporary files.  Run it in its own directory,
    # so multiple instances of mfold won't overwrite each other's
    # files.
    cwd = os.getcwd()
    os.chdir(tempdir)
    
    #seq_file = os.path.join(tempdir, "test.fa")
    seq_file = "test.fa"
    genomefns.write_fasta("sequence", sequence, handle=open(seq_file, 'w'))

    r = mfold(seq_file, NA=NA, T=T)
    output = r.read()
    #print output

    dG = 0
    if output.find("No foldings") < 0:
        ct_file = "%s.ct" % seq_file
        assert os.path.exists(ct_file), "Missing: %s" % ct_file

        # Format:
        #     20 dG =      0.59 sequence
        line = open(ct_file).readline()
        dG = line.strip().split()[3]
        dG = float(dG)

    # Change back to the original directory.
    os.chdir(cwd)

    assert os.path.split(tempdir)[1].startswith("tmp"), tempdir
    os.system("rm -rf %s" % tempdir)

    return dG

def mfold(seq_file, NA=None, NC_CONC=None, T=None):
    # NA        RNA (default) or DNA
    # NC_CONC   Na+ molar concentration (default 1.0)
    # T         Temperature, default 37.
    import pexpect

    cmd = [
        "mfold",
        "SEQ='%s'" % seq_file,
        ]
    if NA is not None:
        assert NA.upper() in ["DNA", "RNA"]
        cmd.append("NA=%s" % NA)
    if NC_CONC is not None:
        cmd.append("NC_CONC=%s" % NC_CONC)
    if T is not None:
        cmd.append("T=%s" % T)

    # mfold writes output directly to TTY, instead of STDOUT.  Need to
    # use pexpect to capture it.
    cmd = " ".join(cmd)
    r = pexpect.spawn(cmd)
    #w, r = os.popen4(cmd)
    #w.close()
    return r
