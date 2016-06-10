import os
import shutil
from pyne import nucname
dirname = "JENDL-AN-2005-linux-patched"
zaiddirname = dirname + "-ZAID/"
if not os.path.exists(zaiddirname):
    os.makedirs(zaiddirname)
files = [ f for f in os.listdir(dirname) if (os.path.isfile(os.path.join(dirname,f)) & f.endswith(".dat")) ]
for file in files:
    nucid = nucname.id(os.path.splitext(file)[0])
    zaid = nucname.zzzaaa(nucid)
    zaidfile = str(zaid) + ".dat"
    shutil.copyfile(os.path.join(dirname , file), os.path.join(zaiddirname, zaidfile ))
