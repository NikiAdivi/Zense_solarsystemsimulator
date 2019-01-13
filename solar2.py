from visual import *
from numpy import linspace,abs,array
import time 
Sun = sphere(radius=695500*40,material=materials.emissive,color=color.orange)
universe=sphere(radius=30066790000,material=materials.texture(data=materials.loadTGA("galaxy_zense")))

mercury = {'a': 57909050.,
           'e': 0.205630,
           'inclination': 7 * pi / 180.,
           'ascension': 0.8436854966,
           'mean_anomaly': 3.0507657193,
           'radius': 2439.7,
           'tilt': 0.1 * pi / 180.,
           'spin': 2 * pi / 4222.6,
           'tex':  materials.texture(data=materials.loadTGA("mercury")),
           'name': "Mercury"}


     
venus = {'a': 108208000.,
         'e': 0.0067,
         'inclination': 3.39 * pi / 180.,
         'ascension': 1.3381895772,
         'mean_anomaly': 0.8746717546,
         'radius': 6051.8,
         'tilt': 177 * pi / 180.,
         'spin': -2 * pi / 2802.,
         'tex':materials.texture(data=materials.loadTGA("venus"),mapping='spherical'),
         'name': "Venus"}

earth = {'a': 149598261.,
         'e': 0.01671123,
         'inclination': 0.,
         'ascension': 0.,
         'mean_anomaly': 6.2398515744,
         'radius': 6378.,
         'tilt': 23 * pi / 180.,
         'spin': 2 * pi / 24.,
         'tex':materials.earth,
         'name': "Earth"}

mars = {'a': 227939100.,
        'e': 0.093315,
        'inclination': 1.85 * pi / 180.,
        'ascension': 0.8676591934,
        'mean_anomaly': 0.3378329113,
        'radius': 3393.5,
        'tilt': 25 * pi / 180.,
        'spin': 2 * pi / 24.66,
        'tex':materials.texture(data=materials.loadTGA("mars"),mapping='spherical'),
        'name': "Mars"}

jupiter = {'a': 778547200.,
           'e': 0.048775,
           'inclination': 1.31 * pi / 180.,
           'ascension': 1.7504400393,
           'mean_anomaly': 0.3284360586,
           'radius': 71400.,
           'tilt': 3 * pi / 180.,
           'spin': 2 * pi / 9.93,
           'tex':materials.texture(data=materials.loadTGA("jupiter"),mapping='spherical'),
           'name': "Jupiter"}

saturn = {'a': 1433449370.,
          'e': 0.055723219,
          'inclination': 2.49 * pi / 180.,
          'ascension': 1.98,
          'mean_anomaly': 5.5911055356,
          'radius': 60000.,
          'tilt': 27 * pi / 180.,
          'spin': 2 * pi / 10.66,
          'tex':materials.texture(data=materials.loadTGA("saturn"),mapping='spherical'),
          'name': "Saturn"}

uranus = {'a': 2876679082.,
          'e': 0.044405586,
          'inclination': 0.77 * pi / 180.,
          'ascension': 1.2908891856,
          'mean_anomaly': 2.4950479462,
          'radius': 25600.,
          'tilt': 98 * pi / 180.,
          'spin': -2 * pi / 17.24,
          'tex':materials.texture(data=materials.loadTGA("uranus"),mapping='spherical'),
          'name': "Uranus"}

neptune = {'a': 4503443661.,
           'e': 0.011214269,
           'inclination': 1.77 * pi / 180.,
           'ascension': 2.3001058656,
           'mean_anomaly': 4.6734206826,
           'radius': 24300.,
           'tilt': 30 * pi / 180.,
           'spin': 2 * pi / 16.11,
           'tex':materials.texture(data=materials.loadTGA("neptune"),mapping='spherical'),
           'name': "Neptune"}

list_of_planets=[mercury,venus,earth,mars,jupiter,saturn,uranus,neptune]

def orbit(m0, e, a, inclination, ascension):
    m = linspace(m0, 2 * pi + m0, 10000)
    ecc_anom = m
    ecc_anom_old = 0
    while 1.e-2 < abs(ecc_anom - ecc_anom_old).max():  # Newton-Raphson solver for eccentric anomaly (E)
        ecc_anom_old = ecc_anom  # assigning previous value for accuracy comparison
        ecc_anom -= (ecc_anom - m - e * sin(ecc_anom)) / (1. - e * cos(ecc_anom))
        # function for E divided by its derivative w.r.t E

    theta = 2. * arctan2(sqrt(1. + e) * sin(ecc_anom / 2.),
                         sqrt(1. - e) * cos(ecc_anom / 2.))  # true anomaly

    r = a * (1 - e * cos(ecc_anom))  # radius

    theasc = theta - ascension

    # conversion to cartesian coordinates:

    x = r * (cos(ascension) * cos(theasc) - sin(ascension) * sin(theasc) * cos(inclination))

    z = r * (sin(ascension) * cos(theasc) + cos(ascension) * sin(theasc) * cos(inclination))

    y = r * (sin(theta - ascension) * sin(inclination))

    coordinates = array((x, y, z))
    return coordinates

def initialise(oneplanet):
    
    planet_obj=sphere(radius=oneplanet['radius']*2000,trail=curve(color=color.white),material=oneplanet['tex']) #2000 is for scaling up 
    planet_obj.rotate(angle=oneplanet['tilt'],axis=vector(0,0,1))
    coordinates=orbit(oneplanet['mean_anomaly'],oneplanet['e'],oneplanet['a'],oneplanet['inclination'],oneplanet['ascension'])
    rotation_rate=365.25*86400*oneplanet['spin']
    timeperiod=(oneplanet['a']/earth['a'])**1.5
    return [planet_obj,coordinates,timeperiod,oneplanet['a'],oneplanet['tilt'],rotation_rate]     #simulation parameters

simulation_planets=[]   #list of simulation parameters of all the planets.

for i in list_of_planets:
    simulation_planets.append(initialise(i))

t=0
dt=1


def update(simulation_oneplanet,t,dt):
    simulation_oneplanet[0].pos= simulation_oneplanet[1][:,(int(-t / simulation_oneplanet[2])) % 10000]
    print simulation_oneplanet[0].pos
    simulation_oneplanet[0].trail.append(simulation_oneplanet[0].pos)#,retain=int(1.5e-9*simulation_oneplanet[3]**1.55/dt))
    simulation_oneplanet[0].rotate(angle=simulation_oneplanet[5]*dt,axis=(-sin(simulation_oneplanet[4]),cos(simulation_oneplanet[4]),0))
    return 0
   

while True: 
    rate(40)
    t=t+100*dt
    for i in simulation_planets:
    	update(i,t,dt)

