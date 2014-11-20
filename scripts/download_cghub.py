#download_cghub.py
import argparse
import os
import subprocess
from genomicode import config

def main():
    parser = argparse.ArgumentParser(description='download_data_from_cghub')
    parser.add_argument(
        '--analysis_id', dest='analysis_id', type=str,
        default=None, help='Given the analysis_id to download')
    parser.add_argument(
        '--analysis_id_file', dest='analysis_id_file', type=str,
        default=None, help='Given the file of analysis_ids to download')
    parser.add_argument(
        '--output', dest='output', default='./', type=str,
        help='path to store the downloaded file')
   
    args = parser.parse_args()
    if args.output:
        assert os.path.exists(args.output),'the path to store the file does not exists'
    assert args.analysis_id or args.analysis_id_file,'please give analysis id or file of analysis ids to download'
    assert not (args.analysis_id and args.analysis_id_file) ,'please give only analysis id or file of analysisids to download'
    key_path = config.cghubkey
    assert os.path.exists(key_path),'cannot find the key file for gtdownload'
    if args.analysis_id_file:
        assert os.path.exists(args.analysis_id_file),'the file %s does not exists'%args.analysis_id_file
    if args.analysis_id_file:
        f=file(args.analysis_id_file)
        text = f.readlines()
        f.close()
        analysis_ids = [i.strip() for i in text]
    elif args.analysis_id:
        analysis_ids = [args.analysis_id]
    for analysis_id in analysis_ids:
        command = ['gtdownload', '-c',key_path, '-d',analysis_id]
        if args.output:
            command.extend(['-p',args.output])
        process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        process.wait()
        error_message = process.communicate()[1]
        if error_message:
            print error_message
        assert os.path.exists(os.path.join(args.output,analysis_id))
    print 'Download finished'



if __name__ == '__main__':
    main()


