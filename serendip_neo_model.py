#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import h5py 
import scipy.interpolate as sio

h=h5py.File("snr_10_diam_e3d.h5","r")
a=n.copy(h["diam_m"].value)
r=n.copy(h["range_m"].value)

rfun = sio.interp1d(a,r)
LD=384400e3


def plot_max_range(D=n.linspace(0.02,100.0,num=100000)):
    plt.loglog(D,rfun(D)/LD)
    plt.xlabel("Diameter (m)")
    plt.ylabel("Maximum range of detection (LD)")
    plt.show()

plot_max_range()


def neo_cumulative_flux(d):
    #(10**(1.568-2.7*n.log10(d)))
    return(36.9828*d**(-2.7))

def neo_density(d):
    return(99.8536*d**(-3.7))

def plot_neo_flux():
    d=10**n.linspace(-2,2,num=100)
    plt.loglog(d,neo_diam_flux(d))
    plt.xlabel("Diameter (m)")
    plt.ylabel("Cumulative number of objects colliding with Earth per year")
    plt.tight_layout()
    plt.show()

def e3d_count():
    d=n.linspace(0.02,10,num=1000000)
    delta_d=n.diff(d)[0]
    rho=neo_density(d)
    R=rfun(d)
    
    A_e3d = 0.5*R**2.0*n.pi/180.0
    A_fb=n.pi*6470e3**2.0



    N_e3d=rho*A_e3d/A_fb

    N_cum=n.sum(delta_d*N_e3d)
    print(N_cum)

    plt.loglog(d,N_e3d)
    plt.show()
    
    plt.loglog(d[0:len(d)],n.cumsum(N_e3d[::-1]*delta_d)[::-1],label="Fixed beam pointing")
    plt.loglog(d[0:len(d)],20*n.cumsum(N_e3d[::-1]*delta_d)[::-1],label="20-position fence scan")
    plt.axvline(1.3/(n.pi*n.sqrt(3.0)),color="black",linestyle="--")
    plt.legend()
    plt.xlabel("Diameter (m)")
    plt.ylabel("Cumulative number of radar detections per year")
    plt.ylim([0.1,2000])
    plt.xlim([0.01,4.8])
    plt.tight_layout()
    plt.savefig("cum_rad_det.png")
#    plt.xlabel(
    plt.show()
    
    plt.loglog(d,N_e3d)
    plt.xlabel("Diameter (m)")
            
    plt.show()
    
#plot_neo_flux()    
e3d_count()

    

