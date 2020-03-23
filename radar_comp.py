#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c

ao_gain=4*n.pi*0.38*n.pi*(305/2.0)**2.0/0.126**2.0
ao_lam=3e8/2.38e9
ao_diam=n.sqrt(ao_gain*ao_lam**2.0/(n.pi**2.0*0.38))
print(ao_diam)
ao={"lam":ao_lam,"D":n.sqrt(ao_gain*ao_lam**2.0/(n.pi**2.0*0.38)),"G":ao_gain,"P":0.9e6,"Tsys":23.0}
dss14={"lam":3e8/8.56e9,"D":70.0,"G":4*n.pi*0.64*n.pi*(70.0/2.0)**2.0/0.035**2.0,"P":0.45e6,"Tsys":18.0}
e3d={"lam":1.3,"G":10**4.3,"D":80.0,"P":5e6*0.25,"Tsys":150.0}
euhf={"lam":0.32,"G":10**4.8,"D":32.0,"P":1.8e6*0.125,"Tsys":70.0}


def comp(r0,r1):
    g0=r0["P"]*r0["G"]**2.0*r0["lam"]**(5/2.0)/r0["Tsys"]
    g1=r1["P"]*r1["G"]**2.0*r1["lam"]**(5/2.0)/r1["Tsys"]

    return(g0/g1)

def range_area(r,D=n.linspace(0.26,10,num=100),B=5.0):
    # beam opening angle
    alpha=r["lam"]/r["D"]
    R10=((r["P"]*r["G"]**2.0*r["lam"]**2.0*0.1*(1.0/4.0)*n.pi*D**2.0)/(10*(4*n.pi)**3.0*c.k*r["Tsys"]*B))**(1/4.0)
    A_cs=0.5*R10**2.0*alpha
    return(A_cs,D)

def neo_flux(D):
    """ flux per m^2 """
    A_e=n.pi*6480e3**2.0
    return(99.854*D**(-3.7)/A_e)

D=n.linspace(0.26,10,num=100)    
A_ao,d=range_area(ao,D)
A_e3d,d=range_area(e3d,D)
A_uhf,d=range_area(euhf,D)
A_dss14,d=range_area(dss14,D)
flux=neo_flux(D)

plt.loglog(d,A_ao,label="Arecibo")

A_ao1 = D*(n.pi/32.0)*n.sqrt(0.1/(c.k*50.0))*n.sqrt(0.9e6*1.0/23.0)*0.38*305.122

plt.loglog(d,A_ao1,label="Arecibo 1")
#plt.loglog(d,A_dss14,label="Goldstone DSS14")
#plt.loglog(d,A_e3d,label="E3D")
#plt.loglog(d,A_uhf,label="UHF")
plt.ylabel("Beam search cross-section area")
plt.xlabel("Object diameter ($m$)")
plt.title("Search cross-section areas")
plt.legend()
plt.show()


plt.loglog(d,flux*A_ao,label="Arecibo")
plt.loglog(d,flux*A_dss14,label="Goldstone DSS14")
plt.loglog(d,flux*A_e3d,label="E3D")
plt.loglog(d,flux*A_uhf,label="UHF")
plt.ylabel("Detections per year / m")
plt.xlabel("Object diameter ($m$)")
plt.title("Detected object flux density")
plt.legend()
plt.show()



print(comp(ao,ao))
print(comp(ao,dss14))
print(comp(ao,e3d))
print(comp(ao,euhf))


print(comp(dss14,ao))
print(comp(dss14,dss14))
print(comp(dss14,e3d))
print(comp(dss14,euhf))


#print(comp(e3d,euhf))
    


