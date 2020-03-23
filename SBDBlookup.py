import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import datetime
import stuffr
import re
import astropy.units as u
from astroquery.jplsbdb import SBDB

import neo_cat

def database_lookup(object_name):
    neo = SBDB.query(object_name)
    elements = neo['orbit']['elements']
    #vecto = [elements['e'], 
     #       elements['a'], 
      #      elements['q'], 
       #     elements['i'],
        #    elements['om'], 
         #   elements['w'], 
          #  elements['ma'], 
           # elements['tp'], 
            #elements['per'], 
    #        elements['n'], 
     #       elements['ad']]
    return(elements)
    # return(vecto)


if __name__ == "__main__":
    neos=neo_cat.read_neos(fname="cneos_closeapproach_data_past.csv")
    names = []
    for n_e_o in neos:
        names.append(n_e_o["name"])
    element_vector = []
    for name in names:

        if name != "66391 Moshup (1999 KW4)":
            element_vector.append(database_lookup(name.replace("(","").replace(")","")))
        if name == "66391 Moshup (1999 KW4)":
            name = "66391 Moshup"
            element_vector.append(database_lookup(name.replace("(","").replace(")","")))
    
    eccentricity = []
    semi_axis = []
    inclination = []
    peri_arg = []

    print("checkpoint")

    for element in element_vector:
        eccentricity.append(float(element['e']))
        semi_axis.append(element['a'].to(u.m)/149597870700.0)
        inclination.append(element['i'].to(u.rad))
        peri_arg.append(element['w'].to(u.rad))
    print(eccentricity[1], type(eccentricity[1]))

    eccentricity = n.array(eccentricity)
    semi_axis = n.array(semi_axis)
    inclination = n.array(inclination)
    peri_arg = n.array(peri_arg)

    plt.rcParams.update({'font.size':16})
    plt.subplot(231)
    plt.plot(semi_axis, eccentricity, '.')
    plt.xlabel("Semi-major axis (AU)")
    plt.ylabel("Eccentricity")
    plt.grid()

    plt.subplot(232)
    plt.plot(semi_axis, inclination/n.pi*180.0, '.')
    plt.xlabel("Semi-major axis (AU)")
    plt.ylabel("Inclination (deg)")
    plt.title("Orbital elements of observable NEO close encounters")
    plt.grid()

    plt.subplot(233)
    plt.plot(semi_axis, peri_arg/n.pi*180.0, '.')
    plt.xlabel("Semi-major axis (AU)")
    plt.ylabel("Argument of perihelion (deg)")
    plt.grid()
    
    plt.subplot(234)
    plt.plot(eccentricity, inclination/n.pi*180.0, '.')
    plt.xlabel("Eccentricity")
    plt.ylabel("Inclination (deg)")
    plt.grid()

    plt.subplot(235)
    plt.plot(eccentricity, peri_arg/n.pi*180.0, '.')
    plt.xlabel("Eccentricity")
    plt.ylabel("Argument of perihelion (deg)")
    plt.grid()
    
    plt.subplot(236)
    plt.plot(inclination/n.pi*180.0, peri_arg/n.pi*180.0, '.')
    plt.xlabel("Inclination (deg)")
    plt.ylabel("Argument of perihelion (deg)")
    plt.grid()
    plt.show()

    print(n.sum(eccentricity)/eccentricity.size)
    print(n.sum(inclination)/inclination.size)
    print(n.sum(semi_axis)/semi_axis.size)
    print(n.sum(peri_arg)/peri_arg.size)

#0.4163445267489712
#0.13019004877165374
#1.580401646090535
#3.3150207207308897



    dates = []
    for objec in neos:
        dates.append(objec["ymd"])
    print(dates[1])

    months = []

    #number per month?
    for dat in dates:
        months.append(datetime.datetime.strptime(dat,'%Y-%b-%d').month)
    
    terrible_implementation = []
    for i in range(1,13):
        terrible_implementation.append(months.count(i))
    plt.plot(range(1,13),terrible_implementation,'.')
    plt.show()