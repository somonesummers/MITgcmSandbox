import sys
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
from glaciome1D import constants, glaciome
import glob
import pickle
import argparse

parser = argparse.ArgumentParser(description='Print the pickel file specified')
parser.add_argument('-f','--file', nargs=1, type=str ,default = None,
                    help='file name to print')
args = parser.parse_args()

files = sorted(glob.glob(args.file[0]))

with open(files[-1], 'rb') as file:
    data = pickle.load(file)
    file.close()

print(data)
