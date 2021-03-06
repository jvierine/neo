#!/usr/bin/env python

import numpy as n
import matplotlib.pyplot as plt
import scipy.constants as c
import datetime
import stuffr
import re
import astropy.units as u

from SBDBlookup import database_lookup
import neo_cat
import neo_horizons

class radar:
    def __init__(self,gain,tx_pwr,duty_cycle,wavelength,noise_temp,max_coh_int_time=0.2):
        self.gain=gain
        self.tx_pwr=tx_pwr
        self.duty_cycle=duty_cycle
        self.wavelength=wavelength
        self.noise_temp=noise_temp
        self.max_coh_int_time=max_coh_int_time
        
class space_object:
    def __init__(self,diameter_m,range_m,spin_period_s=3600.0,radar_albedo=0.1):
        self.diameter_m=diameter_m
        self.range_m=range_m
        self.spin_period_s=spin_period_s
        self.radar_albedo=radar_albedo

def hard_target_s_n(gain_tx, gain_rx,
                    wavelength_m, power_tx,
                    range_tx_m, range_rx_m,
                    diameter_m=0.01,
                    radar_albedo=0.1):
    '''
    Deterine the energy-to-noise ratio for a hard target (signal-to-noise ratio).
    Assume a smooth transition between Rayleigh and optical scattering.
    
    gain_tx - transmit antenna gain, linear
    gain_rx - receiver antenna gain, linear
    wavelength_m - radar wavelength (meters)
    power_tx - transmit power (W)
    range_tx_m - range from transmitter to target (meters)
    range_rx_m - range from target to receiver (meters)
    diameter_m - object diameter (meters)
    radar_albedo - how much is the rcs compared to the rcs of a perfectly 
                 conducting sphere. 
    '''
    
    ##
    ## Determine returned signal power, given diameter of sphere
    ## Ignore resonant regime and use either optical or rayleigh scatter.
    ## This is essentially the same as the NASA SEM model for space
    ## debris RCS. 
    ##
    is_rayleigh = diameter_m < wavelength_m/(n.pi*n.sqrt(3.0))
    is_optical = diameter_m >= wavelength_m/(n.pi*n.sqrt(3.0))
    rayleigh_power = (9.0*power_tx*(((gain_tx*gain_rx)*(n.pi**2.0)*(diameter_m**6.0))/(256.0*(wavelength_m**2.0)*(range_rx_m**2.0*range_tx_m**2.0))))
    optical_power = (power_tx*(((gain_tx*gain_rx)*(wavelength_m**2.0)*(diameter_m**2.0)))/(256.0*(n.pi**2)*(range_rx_m**2.0*range_tx_m**2.0)))
    return(   ((is_rayleigh)*rayleigh_power + (is_optical)*optical_power)*radar_albedo ) 

def incoh_snr_calc(p_s, p_n, epsilon=0.05, B=10.0, t_incoh=3600.0):
    """
    """
    snr=p_s/p_n
    t_epsilon=((p_s+p_n)**2.0)/(epsilon**2.0 * p_s**2.0 * B)

    K=t_incoh*B
#    print(K)
    delta_pn=p_n/n.sqrt(K)
    snr_incoh = p_s/delta_pn
    
    return( snr, snr_incoh, t_epsilon )


def detectability(r,o, t_obs=3600.0,debug=False):
    """
    r = radar
    o = object
    t_obs = observation duration

    returns:
    snr - signal to noise ratio using coherent integration, when doing object discovery with a 
          limited coherent integration duration and no incoherent integration
    snr_incoh - the signal to noise ratio using incoherent integration, when using a priori
                orbital elements to assist in coherent integration and incoherent integration.
                coherent integration length is determined by t_obs (seconds)
    """
    
    doppler_bandwidth=4*n.pi*o.diameter_m/(r.wavelength*o.spin_period_s)

    # for serendipitous discovery
    detection_bandwidth=n.max([doppler_bandwidth,1.0/r.max_coh_int_time, 1.0/t_obs])

    # for detection with a periori know orbit
    # the bandwidth cannot be smaller than permitted by the observation duration.
    incoh_int_bandwidth=n.max([doppler_bandwidth, (1.0/t_obs) ])

    # effective noise power when using just coherent integration 
    p_n0 = c.k*r.noise_temp*detection_bandwidth/r.duty_cycle

    # effective noise power when doing incoherent integration and using a good a priori orbital elements
    p_n1 = c.k*r.noise_temp*incoh_int_bandwidth/r.duty_cycle
    
    p_s=hard_target_s_n(r.gain, r.gain,
                        r.wavelength, r.tx_pwr,
                        o.range_m, o.range_m,
                        diameter_m=o.diameter_m, 
                        radar_albedo=o.radar_albedo)

    snr,snr_incoh,te=incoh_snr_calc(p_s,p_n1,B=incoh_int_bandwidth,t_incoh=t_obs)
    
    snr_coh=p_s/p_n0

    if debug:
        print("discover_snr_coh_%1.2f_s %1.2f track_snr_coh_%1.2f_s %1.2f snr_incoh_%1.1f_s %1.2f required_measurement_time %1.2f dop_bw %1.2f (Hz)"%(1.0/detection_bandwidth,
                                                                                                                                                      snr_coh,
                                                                                                                                                      1.0/incoh_int_bandwidth,
                                                                                                                                                      snr,
                                                                                                                                                      t_obs,
                                                                                                                                                      snr_incoh,
                                                                                                                                                      te,
                                                                                                                                                      doppler_bandwidth))
        
    return(snr_coh,snr_incoh)


def check_radars():
    
    e3d=radar(gain=10**4.3,
              tx_pwr=5e6,
              duty_cycle=0.25,
              wavelength=1.3,
              noise_temp=150.0)
    
    uhf=radar(gain=10**4.8,
              tx_pwr=1.7e6,
              duty_cycle=0.125,
              wavelength=0.32,
              noise_temp=80.0)
    
    arecibo=radar(gain=10**7.5,
                  tx_pwr=1e6,
                  duty_cycle=1.0,
                  wavelength=0.125,
                  noise_temp=50.0)
    
    dists=10**n.linspace(2,9,num=100)
    
    o = space_object(diameter_m=2,
                     range_m=3.6e5*1e3,
                     spin_period_s=5*60,
                     radar_albedo=0.1)

    print("EISCAT")
    snr_coh,snr_incoh=detectability(uhf,o, t_obs=3600.0)
    print(snr_coh)
    print(snr_incoh)

    print("Arecibo")
    snr_coh,snr_incoh=detectability(arecibo,o, t_obs=3600.0)
    print(snr_coh)
    print(snr_incoh)


def arecibo_radar():
    return(radar(gain=10**7.5,
                 tx_pwr=1e6,
                 duty_cycle=1.0,
                 wavelength=0.125,
                 noise_temp=50.0))

def e3d_radar():
    return(radar(gain=10**4.3,
                 tx_pwr=5e6,
                 duty_cycle=0.25,
                 wavelength=1.3,
                 noise_temp=150.0))
def uhf_radar():
    return(radar(gain=10**4.8,
                 tx_pwr=1.8e6,
                 duty_cycle=0.125,
                 wavelength=0.32,
                 noise_temp=90.0))
    

    
def catalog_check(r,
                  planar_array=True,
                  time_window=2*24*3600.0,
                  fname="cneos_closeapproach_data_past.csv",
                  debug=False):

    """ 
    go through catalog
    check is object can be detected at closest distance
    if yes, check using horizons if it is above horizon and at a detectable range.
    """
    gain0=r.gain
    
    usable_list = []
    counter = 0
    LD=384400e3
    neos=neo_cat.read_neos(fname=fname)
    det_r=[]
    det_d=[]
    obs_time_times = []

    n_detectable=0
    for neo in neos:

        o = space_object(diameter_m=0.5*(neo["d_min"]+neo["d_max"]),
                         range_m=LD*neo["dist_ld"],
                         spin_period_s=60*60*0.005*0.5*(neo["d_min"]+neo["d_max"]),
                         radar_albedo=0.1)
        
        snr_coh,snr_incoh=detectability(r,o, t_obs=3600.0)

        name=re.search(".*\((.*)\)",neo["name"]).group(1)

        serendipitously_discoverable = []


        
        if snr_incoh > 10.0:
            # we may have a chance. let's look at elevation and observability using a real ephemeris
            if False:
                print("# %s min_dist_ld %1.2f max_snr_coh %1.2f max_snr_incoh %1.2f"%(name,neo["dist_ld"],snr_coh,snr_incoh))

            # I hate datetime calculations. I'm pretty sure this only works if your computer timezone is UTC
            # this is not production code, so I'm not going to bother to figure out how to do this more
            # generally.
            d0=datetime.datetime.strptime("%s %s"%(neo["ymd"],neo["hour"]), "%Y-%b-%d %H:%M")
            timestamp = (d0-datetime.datetime(1970,1,1)).total_seconds()
            d0=stuffr.unix2date(timestamp-time_window)
            d1=stuffr.unix2date(timestamp+time_window)
            start_t=d0.strftime('%Y-%m-%d %H:%M')
            stop_t=d1.strftime('%Y-%m-%d %H:%M')

            # get ephemeris during an observing window around closest approach
            h_ranges,h_range_rates,h_els,h_dates= neo_horizons.check_detectability(obj_id=name,
                                                                                   start=start_t,
                                                                                   stop=stop_t,
                                                                                   step="5m")
            obs_time = []
            n_obs=len(h_ranges)
            coh_max = 0
            max_snr=0
            min_dist=0
            max_el=0
            max_date=""
            total_obs_time=0.0
            total_s_obs_time=0.0            
            for oi in range(n_obs):
                o.range_m=h_ranges[oi]
                if planar_array:
                    r.gain=gain0*n.sin(n.pi*h_els[oi]/180.0)
                
                m_snr_coh,m_snr_incoh=detectability(r,o, t_obs=3600.0)

                if snr_incoh > 10:
                    total_obs_time+=5.0
                    if m_snr_incoh > max_snr:
                        max_snr=m_snr_incoh
                        min_dist=h_ranges[oi]/LD
                        max_el=h_els[oi]
                        max_date=h_dates[oi]
                        #                    print("%s dist_ld %1.2f snr_coh %1.2f snr_incoh/hour %1.2f"%(name,h_ranges[oi]/LD,m_snr_coh,m_snr_incoh))
                    obs_time.append(h_dates[oi])
                if m_snr_coh > 10.0:
                    coh_max=m_snr_coh
                    total_s_obs_time+=5.0
            if max_snr > 10.0:
                n_detectable+=1
                print("%d -> %s %s min_dist_ld %1.2f elevation %1.2f snr_coh %1.2f snr_incoh/hour %1.2f obs_time %1.2f (min) coh_obs_time %1.2f (min)"%(n_detectable,name,max_date,min_dist,max_el,coh_max,max_snr,total_obs_time,total_s_obs_time))
                counter +=1
                usable_list.append([name, max_date, min_dist, max_el, max_snr,neo["d_max"]]) #saves the observable ones
                
                obs_time_times.append(obs_time)
            if debug:
                print("----")
        else:
            if debug:
                print("%s not observable"%(name))
    print(counter)

    print(len(obs_time_times), len(obs_time_times[0]))


    object_observability_durations = []
    for time in obs_time_times:
        duration = 0
        for step in time:
            duration +=1

    ranges = []
    diameters = []
    snrs = []
    dates = []
    elements = []
    for i in range(len(usable_list)):
        ranges.append(usable_list[i][2])
        snrs.append(usable_list[i][4])
        diameters.append(usable_list[i][5])
        dates.append(usable_list[i][1])
        elements.append(database_lookup(usable_list[i][0]))
    #plt.plot(ranges, diameters,'.')
    #plt.grid()
    #plt.title("Diameter as a function of range of observable NEOs")
    #plt.xlabel("Range (LD)")
    #plt.ylabel("Diameter")
    #plt.show()

    plt.plot(ranges,snrs,'.')
    plt.show()


    date_attempt = datetime.datetime.strptime(dates[0][:-6],'%Y-%b-%d')
    months = []

    #number per month?
    for dat in dates:
        months.append(datetime.datetime.strptime(dat[:-6],'%Y-%b-%d').month)
    
    terrible_implementation = []
    for i in range(1,13):
        terrible_implementation.append(months.count(i))

    print(months)
    print(terrible_implementation)
    
    print(terrible_implementation)
    date_axis = ['2019-mar', '2019-apr', '2019-may', '2019-jun','2019-jul','2019-aug','2019-sep','2019-oct', '2019-nov','2019-dec','2020-jan', '2020-feb','2020-mar']
    proper_date_axis = []
    for dat in date_axis:
        proper_date_axis.append(datetime.datetime.strptime(dat,'%Y-%b'))
    
    plt.plot(range(1,13), terrible_implementation, '.')
    plt.show()

    eccentricity = []
    semi_axis = []
    inclination = []
    peri_arg = []

    eccentricity2 = []
    semi_axis2 = []
    inclination2 = []
    peri_arg2 = []

    for element in elements:
        eccentricity.append(float(element['e']))
        semi_axis.append(element['a'].to(u.m)/149597870700.0)
        inclination.append(element['i'].to(u.rad))
        peri_arg.append(element['w'].to(u.rad))

    element_vector = []
    diams = []
    for neo in neos:
        diams.append(0.5*neo["d_min"]+0.5*neo["d_max"])
        name = neo["name"]
        if name != "66391 Moshup (1999 KW4)":
            element_vector.append(database_lookup(name.replace("(","").replace(")","")))
        if name == "66391 Moshup (1999 KW4)":
            name = "66391 Moshup"
            element_vector.append(database_lookup(name.replace("(","").replace(")","")))

    eccentricity = n.array(eccentricity)
    semi_axis = n.array(semi_axis)
    inclination = n.array(inclination)
    peri_arg = n.array(peri_arg)

    for element in element_vector:
        eccentricity2.append(float(element['e']))
        semi_axis2.append(element['a'].to(u.m)/149597870700.0)
        inclination2.append(element['i'].to(u.rad))
        peri_arg2.append(element['w'].to(u.rad))

    eccentricity2 = n.array(eccentricity2)
    semi_axis2 = n.array(semi_axis2)
    inclination2 = n.array(inclination2)
    peri_arg2 = n.array(peri_arg2)

    plt.rcParams.update({'font.size':16})
    plt.rcParams["figure.figsize"] = (30,10)
    plt.subplot(231)
    plt.plot(semi_axis2, eccentricity2, '.')
    plt.plot(semi_axis, eccentricity, 'x', color='r', markersize = 12)
    plt.xlabel("Semi-major axis (AU)")
    plt.ylabel("Eccentricity")
    plt.grid()

    plt.subplot(232)
    plt.plot(semi_axis2, inclination2/n.pi*180.0, '.')
    plt.plot(semi_axis, inclination/n.pi*180.0, 'x', color = 'r', markersize = 12)
    plt.xlabel("Semi-major axis (AU)")
    plt.ylabel("Inclination (deg)")
    plt.title("Orbital elements of observable NEO close encounters and total population, EISCAT 3D")
    plt.grid()

    plt.subplot(233)
    plt.plot(semi_axis2, peri_arg2/n.pi*180.0, '.')
    plt.plot(semi_axis, peri_arg/n.pi*180.0, 'x',color='r', markersize = 12)
    plt.xlabel("Semi-major axis (AU)")
    plt.ylabel("Argument of perihelion (deg)")
    plt.grid()
    
    plt.subplot(234)
    plt.plot(eccentricity2, inclination2/n.pi*180.0, '.')
    plt.plot(eccentricity, inclination/n.pi*180.0, 'x', color='r', markersize = 12)
    plt.xlabel("Eccentricity")
    plt.ylabel("Inclination (deg)")
    plt.grid()

    plt.subplot(235)
    plt.plot(eccentricity2, peri_arg2/n.pi*180.0, '.')
    plt.plot(eccentricity, peri_arg/n.pi*180.0, 'x', color='r', markersize = 12)
    plt.xlabel("Eccentricity")
    plt.ylabel("Argument of perihelion (deg)")
    plt.grid()
    
    plt.subplot(236)
    plt.plot(inclination2/n.pi*180.0, peri_arg2/n.pi*180.0, '.')
    plt.plot(inclination/n.pi*180.0, peri_arg/n.pi*180.0, 'x',color='r', markersize = 12)
    plt.xlabel("Inclination (deg)")
    plt.ylabel("Argument of perihelion (deg)")
    plt.grid()
    plt.tight_layout()
    plt.show()

    print(len(inclination))




    


#                    print("%s dist_ld %1.2f snr_coh %1.2f snr_incoh %1.2f diam %1.1f-%1.1f m"%(neo["name"],neo["dist_ld"],snr_coh,snr_incoh,neo["d_min"],neo["d_max"]))
#det_r.append(neo["dist_ld"]*LD)
 #               det_d.append(0.5*(neo["d_min"]+neo["d_max"]))
  #  det_r=n.array(det_r)
   # det_d=n.array(det_d)
    
#    plt.loglog(det_r/1e3,det_d,".")
#    plt.xlabel("Distance (km)")
#    plt.ylabel("Diameter (m)")
#    plt.title("Feasible NEO tracks within MPC catalog")
#    plt.show()
        
        

if __name__ == "__main__":
    e3d=e3d_radar()
    uhf=uhf_radar()    

    print("Past year")
    #print("EISCAT UHF")
    #catalog_check(uhf,planar_array=False,fname="cneos_closeapproach_data_past.csv")
    print("EISCAT 3D")
    catalog_check(e3d,planar_array=True,fname="cneos_closeapproach_data_past.csv")
    
    #print("One year into future")
    #print("EISCAT UHF")
    #catalog_check(uhf,planar_array=False,fname="cneos_closeapproach_data_future.csv")
    #print("EISCAT 3D")
    #catalog_check(e3d,planar_array=True,fname="cneos_closeapproach_data_future.csv")
    
