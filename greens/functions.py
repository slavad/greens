import greens.raw_functions as grf
from greens.raw_functions import testgreenfun01, rootslipunova

def greenlip(R,R_,t,R_in,R_out,n,root):
    return grf.greenlip(R,R_,t,R_in,R_out,n,root, root.shape[0])