import argparse 
import os
import shutil
import sys

OPTION = 'option'
H5 = 'h5'
TSV = 'tsv'

def parse_input():
    '''
    Parses the commandline input, returns a dict
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('file_list', metavar='File', nargs='+')
    parser.add_argument('-t', required=True, dest=OPTION)
    args = parser.parse_args()
    return vars(args)


if __name__ == '__main__':
    
    args = parse_input()
    file_list = args['file_list']
    samples = ['.'.join(os.path.basename(x).split('.')[:-1]) for x in file_list]
    for i, f in enumerate(file_list):
        # file paths are like /cromwell_root/.../.../<samplename>.<file type>.<file extension> where file type is 'abundance' or run_info, etc.
        bn = os.path.basename(f)
        contents = bn.split('.')
        if args[OPTION] == H5:
            file_extension = contents[-1]
            file_type = contents[-2]
            sample_name = '.'.join(contents[:-2])
            if not os.path.exists(sample_name):
                os.mkdir(sample_name)
            shutil.move(f, os.path.join(sample_name, '%s.%s' % (file_type, file_extension)))
        elif args[OPTION] == TSV:
            file_extension = '.'.join(contents[-3:])
            sample_name = '.'.join(contents[:-3])
            if not os.path.exists(sample_name):
                os.mkdir(sample_name)
            shutil.move(f, os.path.join(sample_name, file_extension))
        else:
            sys.stdout.write('Please specify a valid option')
            sys.exit(1)
