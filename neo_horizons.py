#!/usr/bin/env python

import numpy as n
import scipy.constants as c
from astroquery.jplhorizons import Horizons
import neo_snr

def check_detectability(obj_id="2020 DA4",
                        obs_lat=69.3908,
                        obs_lon=20.2673,
                        obs_el=0.0,
                        start="2020-01-01",
                        stop="2020-06-01",
                        step="2h",
                        min_el=30.0,
                        id_type="smallbody",
                        debug=False):
    
    e3d = {'lon': obs_lon, 'lat': obs_lat, 'elevation': obs_el}
    obj = Horizons(id=obj_id,
                   location=e3d,
                   epochs={"start":start,
                           "stop":stop,
                           "step":step},id_type=id_type)

    t=obj.ephemerides(quantities="4,20",get_raw_response=False)
    if debug:
        print(t)
    els=[]
    ranges=[]
    range_rates=[]
    dates=[]
    for r in t:

        name=r[0]
        date=r[1]
        az=r[7]   # deg 
        el=r[8]   # deg
        
        if el > min_el:
            delta=r[9]*c.au   # m
            delta_rate=r[10]*1e3  #m/s
            #            print("el %1.2f range %1.2g range-rate %1.2f gain %1.2f"%(el,delta,delta_rate/1e3,g_dB))
            els.append(float(el))
            ranges.append(delta)
            range_rates.append(delta_rate)
            dates.append(r[1])
    return(n.array(ranges),n.array(range_rates),n.array(els),dates)

if __name__ == "__main__":
    check_detectability(obj_id="2002 PZ39",id_type="smallbody",debug=True)
