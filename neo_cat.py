#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt


def read_neos(fname="cneos_closeapproach_data.csv"):

    n_list=[]
    
    f=file(fname,"r")
    for l in f.readlines():
        l=l.split(",")
        name=l[0].strip("\"")
        if name == "Object":
            continue
#        print(name)
        date=l[1].split(" ")
        ymd=date[0].strip("\"")
        hour=date[1][0:5]
 #       print("%s %s %s"%(name,ymd,hour))
#        print(l)
        dist_ld=float(l[2].split(" ")[0].strip("\""))
 #       print(dist_ld)
        diam=l[7]
#        print(diam)
        diam=diam.split(" ")
#        print(diam)
        d_min=float(diam[0].strip("\""))
        n_par=len(diam)
        d_max=float(diam[n_par-2])
#        print("%s %1.2f %1.2f-%1.2f"%(name,dist_ld,d_min,d_max))
        n_list.append({"name":name,"date":date,"d_min":d_min,"d_max":d_max,"dist_ld":dist_ld})
    return(n_list)



if __name__ == "__main__":
    neos=read_neos()
    print(len(neos))
    
