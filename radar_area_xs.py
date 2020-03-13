#!/usr/bin/env python


import numpy as n
import matplotlib.pyplot as plt


def collecting_area_fraction(r=10000.0,R_e=6380.0,theta=n.pi/180.0):
    """
    Estimate what is the collecting area of the E3D radar beam up to height of r km
    compared with the collecting area of Earth (fireball statistics)
    """
    Ap=0.5*r**2.0*theta
    A=n.pi*R_e**2.0
    f=Ap/A
    return(f)


if __name__ == "__main__":

    print(collecting_area_fraction())

    
