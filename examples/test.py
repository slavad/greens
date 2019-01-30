import greens
import pdb
import numpy as np

def test_green_fun_01(R_in,R_out,R_,n,t,n_root,n_R):
    dR=(R_out-R_in)/(n_R-1)
    l=1.0/(4.0-2*n)
    x_out = 1.0
    x_in = (R_in/R_out)**(1.0-n/2.0)
    root = greens.rootslipunova(n_root,l,x_in,x_out)
    R=R_in
    while R <= R_out:
        Sigma,Mdot = greens.greenlip(R,R_,t,R_in,R_out,n,root, root.shape[0])
        print(R,Sigma, Mdot)
        R=R+dR

def main():
    R_in = 1.0
    R_out = 10.0
    R_ = 2.0
    n = 0.75
    t = 1.1
    n_root = 2000
    n_R = 200
    test_green_fun_01(R_in,R_out,R_,n,t,n_root,n_R)
    #greens.testgreenfun01(R_in,R_out,R_,n,t,n_root,n_R)

main()