import numpy as np
import ctypes

rsfind = ctypes.cdll.LoadLibrary("./rsfind.so")

genoFilePath = "/homes/pappleby/data/ukb/geno/chr22.bgen"
dbFilePath = "/homes/pappleby/data/ukb/geno/chr22.bgen.bgi"
rsid = "rs143816252"

rtn = rsfind.init_rsid_search(genoFilePath, dbFilePath, rsid);
print "RTN", rtn
