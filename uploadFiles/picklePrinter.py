import sys
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
sys.path.append('/Users/psummers8/Documents/glaciome1D')
from glaciome1D import constants, glaciome
import glob
import pickle
import argparse

parser = argparse.ArgumentParser(description='Print the pickel file specified')
parser.add_argument('-f','--file', nargs=1, type=str ,default = None,
                    help='file name to print')
parser.add_argument('filename',default = None, type = str)
args = parser.parse_args()

if(args.filename != None):
    fileName = args.filename
else:
    fileName = args.file[0]
files = sorted(glob.glob(fileName))

with open(files[-1], 'rb') as file:
    data = pickle.load(file)
    file.close()

print(data)
