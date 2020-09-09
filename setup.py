import humann
import os

path = humann.__path__

print("Found Humann installation  at:", path)
print("Replace humann scripts with parallel humann")

os.system("cp parallel_humann/search/* {}/search".format(path[0]))
os.system("cp parallel_humann/*.py {}".format(path[0]))

print("Installation complete!")

os.system("humann --version")
