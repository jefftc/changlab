#R.py
import subprocess
import argparse

def run_R(script):
    command=['R','CMD','BATCH']
    command.append(script)
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    process.wait()
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
