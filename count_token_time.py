"""

что нужно ?
найти токены
распарсить время
занести в таблицу

"""
import argparse
import glob

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir")
args = parser.parse_args()
dir = args.dir()
dir = "./"
for token in glob.iglob(dir + "/**/token**/*"):
    print(token)
