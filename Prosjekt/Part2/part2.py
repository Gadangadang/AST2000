import numpy as np
import matplotlib.pyplot as plt

from ast2000tools.space_mission import SpaceMission
import ast2000tools.utils as utils
from ast2000tools.constants import *
from ast2000tools.solar_system import SolarSystem

utils.check_for_newer_version()
# Construct SpaceMission instance for my mission
seed = utils.get_seed('Sgfrette')
mission = SpaceMission(seed)
system = SolarSystem(seed)
"""
def numeric_analytic_orbit():

    e = [i for i in system.eccentricities]
    a = [i for i in system.semi_major_axes]
    print(len(a),len(e))
    name = [i for i in range(8)]
    x= zip(e,a,name)

    def orbit_analytical(e,a,f):
        r = (a*(1-e**2))/(1+e*np.cos(f))
        return r

    angle = np.linspace(0,2*np.pi,1000)

    for i in x:
        r = orbit_analytical(i[0],i[1],angle)
        x = r*np.cos(angle)
        y = r*np.sin(angle)
        plt.plot(x, y, "--",label="{}".format(i[2]))


    Planets = np.zeros((8,5))
    for i in range(len(Planets)):
        Planet0_x0 = [system.initial_positions[0,i],system.initial_positions[1,i]]
        Planet0_v0 = [system.initial_velocities[0][i],system.initial_velocities[1][i]]
        Planets[i,0] = Planet0_x0[0]
        Planets[i,1] = Planet0_x0[1]
        Planets[i,2] = Planet0_v0[0]
        Planets[i,3] = Planet0_v0[1]
        Planets[i,4] = i

    Gr = G_sol
    sunmass = system.star_mass


    def integrator(Planet):
        def gravity(x,y):
            grav = -(Gr*sunmass)/(np.sqrt(x**2+y**2))**2
            ax = grav*x /(np.sqrt(x**2+y**2))
            ay = grav*y/(np.sqrt(x**2+y**2))
            return np.asarray([ax,ay])
        x0x,x0y,v0x,v0y,planet_index = Planet

        N = 1400000
        time = 7.6*40 #time for 20 rotations for our home planet
        dt = time/N
        t = dt
        x = np.zeros((N,2),float)
        x[0,:] = [x0x,x0y]
        v = np.asarray([v0x,v0y])
        a_i = gravity(x[0,0],x[0,1])
        for i in range(N-1):
            t+=dt
            x[i+1,:] = x[i,:]+(v*dt + 0.5*a_i*dt**2)
            a_iplus1 = gravity(x[i+1,0],x[i+1,1])
            v += 0.5*( a_i + a_iplus1 )*dt
            a_i = a_iplus1
        return x, planet_index


    def plot_func(Planet):
        planet_orbit,planet_index = integrator(Planet)
        plt.plot(planet_orbit[:,0], planet_orbit[:,1],label="{}".format(planet_index))

    for Planet in Planets:
        plot_func(Planet)


    plt.legend()
    plt.xlabel("x-position (AU)")
    plt.ylabel("y-position (AU)")
    plt.savefig("Numeric-analy-moresteps.jpeg")

    plt.show()
"""
def star_planets_cm_movements():
    sunmass = system.star_mass


    Gr = G_sol
    def integrator(planets,sun,planets_index):

        def gravity_on_planet(x_sun,y_sun,x_planet,y_planet,index):

            r = np.asarray([x_sun-x_planet, y_sun-y_planet])
            grav = (Gr*sunmass*masses[index])/(np.linalg.norm(r))**2
            ax = grav*r[0] /(np.linalg.norm(r))
            ay = grav*r[1]/(np.linalg.norm(r))

            return np.asarray([ax,ay],float)

        number_of_planets = int(len(planets))
        cm = np.array([0,0])
        masses = [system.masses[i] for i in planets_index]
        x0xsun,x0ysun,v0xsun,v0ysun = sun

        N = 2800000
        time = 7.6*40

        dt = time/N
        t = dt
        x_planets = np.asarray([np.zeros((N,2),float) for i in range(number_of_planets)])
        v_planets = [np.asarray([planets[i,2],planets[i,3]]) for i in range(number_of_planets)]
        for i in range(number_of_planets):
            x_planets[i,0,:] = [planets[i,0],planets[i,1]]
        a_i_planets = np.zeros((number_of_planets,2))
        a_iplus1_planets = np.zeros((number_of_planets,2))

        for i in range(number_of_planets):
            a_i_planets[i] = gravity_on_planet(x0xsun,x0ysun,planets[i,0],planets[i,1],i)

        x_sun = np.zeros((N,2),float)
        x_sun[0,:] = [x0xsun,x0ysun]
        v_sun = np.asarray([v0xsun,v0ysun],float)
        a_i_sun = -sum(a_i_planets)/sunmass


        for i in range(N-1):
            t+=dt
            cm_x = [x_planets[g,i,0]*masses[g] for g in range(number_of_planets)]
            cm_y = [x_planets[g,i,1]*masses[g] for g in range(number_of_planets)]
            cm = [(sum(cm_x)+sunmass*x_sun[i,0])/(sum(masses)+sunmass), \
            (sum(cm_y)+sunmass*x_sun[i,1])/(sum(masses)+sunmass)]

            for j in range(number_of_planets):
                x_planets[j,i+1,:] = x_planets[j,i,:]+(v_planets[j]*dt + 0.5*a_i_planets[j]/masses[j]*dt**2)-cm
            x_sun[i+1,:] = x_sun[i,:]+(v_sun*dt + 0.5*a_i_sun*dt**2)-cm

            for k in range(number_of_planets):
                a_iplus1_planets[k] = gravity_on_planet(x_sun[i+1,0],x_sun[i+1,1],x_planets[k,i+1,0],\
                x_planets[k,i+1,1],k)
            a_iplus1_sun = -sum(a_i_planets)/sunmass

            for l in range(number_of_planets):
                v_planets[l] += 0.5*(a_i_planets[l] + a_iplus1_planets[l])/masses[l]*dt

            a_i_planets = a_iplus1_planets
            v_sun += 0.5*( a_i_sun + a_iplus1_sun )*dt
            a_i_sun = a_iplus1_sun



        return x_planets,x_sun

    """Input index of planets wanted in simulation"""
    planets_index = [0,1,2,3,4,5,6,7]
    planets = []
    for i in planets_index:
        planets.append([system.initial_positions[0,i],system.initial_positions[1,i],\
        system.initial_velocities[0][i],system.initial_velocities[1][i]])
    planets= np.asarray(planets)
    Sun = [0,0,0,0]

    planet_orbit, sun_orbit = integrator(planets,Sun,planets_index)
    for i in range(len(planets_index)):
        plt.plot(planet_orbit[i,:,0], planet_orbit[i,:,1], label="Planet nr.{}".format(planets_index[i]))



    plt.plot(sun_orbit[:,0], sun_orbit[:,1], label="Sun")
    plt.xlabel("x-position (AU)")
    plt.ylabel("y-position (AU)")
    plt.legend(loc='lower right')
    plt.savefig("garvity-sun-pos.jpeg")
    plt.show()

    plt.plot(sun_orbit[:,0], sun_orbit[:,1], label="Sun")
    plt.xlabel("x-position (AU)")
    plt.ylabel("y-position (AU)")
    plt.legend(loc='lower left')
    plt.savefig("sun-all-grav.jpeg")
    plt.show()

star_planets_cm_movements()
"""
def write_radial_vel():

    sunmass = system.star_mass
    planetmass = system.masses[0]
    my = sunmass*planetmass/(sunmass+planetmass)
    masstot = sunmass+planetmass


    sigma = 0.0000126*1/5
    N = 152000
    noise = np.random.normal(0,sigma, size = int(N))
    V_rad = [noise[0]]
    Vofsystem = 46 # AU/year

    Gr = G_sol
    def integrator(planet,sun):

        def gravity_on_sun(x_sun,y_sun,x_planet,y_planet):

            r = np.asarray([x_sun-x_planet, y_sun-y_planet])
            grav = -(Gr*sunmass*planetmass)/(np.linalg.norm(r))**2
            ax = grav*r[0] /(np.linalg.norm(r))
            ay = grav*r[1]/(np.linalg.norm(r))

            return np.asarray([ax,ay],float)

        x0x,x0y,v0x,v0y = planet
        x0xsun,x0ysun,v0xsun,v0ysun = sun

        N = 152000
        time = 7.6*2 #time for 20 rotations for our home planet
        #print(time)
        dt = time/N
        t = dt

        x_planet = np.zeros((N,2),float)
        x_sun = np.zeros((N,2),float)

        cm = np.array([0,0])

        x_planet[0,:] = [x0x,x0y]
        x_sun[0,:] = [x0xsun,x0ysun]

        v_planet = np.asarray([v0x,v0y],float)
        v_sun = np.asarray([v0xsun,v0ysun],float)
        F_sun0 = gravity_on_sun(x0xsun,x0ysun,x0x,x0y)
        a_i_sun = F_sun0/sunmass
        a_i_planet = -F_sun0/planetmass
        for i in range(N-1):
            t+=dt
            cm = np.asarray([(sunmass*x_sun[i,0] + planetmass*\
            x_planet[i,0])/(sunmass+planetmass),\
            (sunmass*x_sun[i,1] + planetmass*\
            x_planet[i,1])/(sunmass+planetmass)])

            x_planet[i+1,:] = x_planet[i,:]+(v_planet*dt + 0.5*a_i_planet*dt**2)-cm
            x_sun[i+1,:] = x_sun[i,:]+(v_sun*dt + 0.5*a_i_sun*dt**2)-cm
            F_sun = gravity_on_sun(x_sun[i+1,0],x_sun[i+1,1],x_planet[i+1,0],\
            x_planet[i+1,1])
            a_iplus1_sun = F_sun/sunmass
            a_iplus1_pl = -F_sun/planetmass
            v_planet += 0.5*( a_i_planet + a_iplus1_pl)*dt
            a_i_planet = a_iplus1_pl
            v_sun += 0.5*( a_i_sun + a_iplus1_sun )*dt
            a_i_sun = a_iplus1_sun

            v_inp = np.linalg.norm(v_sun)*x_sun[i+1,1]/np.linalg.norm(x_sun[i+1])
            V_rad.append(v_inp+noise[i+1])
        return x_planet,x_sun, V_rad, dt

    Planet0 = [system.initial_positions[0,0],system.initial_positions[1,0],\
    system.initial_velocities[0][0],system.initial_velocities[1][0]]

    Sun = [0,0,0,0]

    planet_orbit, sun_orbit, V_rad, dt = integrator(Planet0,Sun)

    V_rad = np.asarray(V_rad) + Vofsystem


    infile = open("Radial_vel.txt","w")
    for i in range(len(V_rad)):
        infile.write("{} {}\n".format(i*dt,V_rad[i]))
    infile.close()

def read_radial_vel():
    L = 10
    sigma = 0.00001/5
    sigma_P = 7.5/5
    sigma_t0 = 5.35/5
    max_val_v = 0.000011+sigma
    min_val_v = 0.000011-sigma
    max_val_P = 7.5+sigma_P
    min_val_P = 7.5-sigma_P
    max_val_t0 = 5.35+sigma_t0
    min_val_t0 = 5.35-sigma_t0
    val_of_cm = 46


    infile = open("Radial_vel.txt", "r")
    infile.readline()

    time = []
    vel = []
    for line in infile:
        word = line.split()
        time.append(float(word[0]))
        vel.append(float(word[1]))


    time = np.asarray(time)
    vel_imp = np.asarray(vel)



    for i in range(len(vel_imp)):
        vel_imp[i] -= val_of_cm
    N = len(time)



    def minste_kvadrat(P,v_star, t_0, vel_imp):
        value = 0
        for i in range(len(time)):
            a = (vel_imp[i]-v_star*np.cos((2*np.pi/P)*(time[i]-t_0)))**2
            value += a
        return P,v_star, t_0,value


    d_value_P = (max_val_P-min_val_P)/(L-1)
    d_value_v = (max_val_v-min_val_v)/(L-1)
    d_value_t0 = (max_val_v-min_val_v)/(L-1)

    P_values = np.asarray([min_val_P + d_value_P*i for i in range(L)])
    v_star_values = np.asarray([min_val_v + d_value_v*i for i in range(L)])
    t_0_values = np.asarray([min_val_t0 + d_value_t0*i for i in range(L)])

    minste_kvadrat_liste = []

    for n in range(len(P_values+1)):
        for g in range(len(P_values+1)):
            for j in range(len(P_values+1)):
                b = minste_kvadrat(P_values[n],v_star_values[g],t_0_values[j],vel_imp)
                minste_kvadrat_liste.append(b)

    minste_kvadrat_liste = np.asarray(minste_kvadrat_liste)

    for i in range(len(minste_kvadrat_liste)):
        if minste_kvadrat_liste[i,3] == min(minste_kvadrat_liste[:,3]):
            index = np.where(minste_kvadrat_liste[i,3])

    #index = np.where(min(minste_kvadrat_liste[:,3]))


    print(minste_kvadrat_liste[index[0]])
    P,v_star,t_0,value =minste_kvadrat_liste[index[0]][0]

    def V(t):
        return v_star*np.cos((2*np.pi*(t-t_0)/P))

    plt.plot(time,vel_imp,label="Data")
    plt.plot(time,V(time),label="Estimated curve")
    plt.plot()

    plt.xlabel("Years")
    plt.ylabel("AU/Years")
    plt.legend()
    plt.savefig("read-vals.jpeg")
    plt.show()

def simulate_eclipse():
    sunmass = system.star_mass
    planetmass = system.masses[0]



    Gr = G_sol
    def integrator(planet,sun):

        def gravity_on_sun(x_sun,y_sun,x_planet,y_planet):

            r = np.asarray([x_sun-x_planet, y_sun-y_planet])
            grav = -(Gr*sunmass*planetmass)/(np.linalg.norm(r))**2
            ax = grav*r[0] /(np.linalg.norm(r))
            ay = grav*r[1]/(np.linalg.norm(r))

            return np.asarray([ax,ay],float)

        x0x,x0y,v0x,v0y = planet
        x0xsun,x0ysun,v0xsun,v0ysun = sun

        N = 76000
        time = 7.6 #time for 20 rotations for our home planet
        #print(time)
        dt = time/N
        t = dt

        x_planet = np.zeros((N,2),float)
        x_sun = np.zeros((N,2),float)

        cm = np.array([0,0])

        x_planet[0,:] = [x0x,x0y]
        x_sun[0,:] = [x0xsun,x0ysun]

        v_planet = np.asarray([v0x,v0y],float)
        v_sun = np.asarray([v0xsun,v0ysun],float)

        F_sun0 = gravity_on_sun(x0xsun,x0ysun,x0x,x0y)
        a_i_sun = F_sun0/sunmass
        a_i_planet = -F_sun0/planetmass

        for i in range(N-1):
            t+=dt
            cm = np.asarray([(sunmass*x_sun[i,0] + planetmass*\
            x_planet[i,0])/(sunmass+planetmass),\
            (sunmass*x_sun[i,1] + planetmass*\
            x_planet[i,1])/(sunmass+planetmass)])

            x_planet[i+1,:] = x_planet[i,:]+(v_planet*dt + 0.5*a_i_planet*dt**2)-cm
            x_sun[i+1,:] = x_sun[i,:]+(v_sun*dt + 0.5*a_i_sun*dt**2)-cm

            F_sun = gravity_on_sun(x_sun[i+1,0],x_sun[i+1,1],x_planet[i+1,0],\
            x_planet[i+1,1])
            a_iplus1_sun = F_sun/sunmass
            a_iplus1_pl = -F_sun/planetmass
            v_planet += 0.5*( a_i_planet + a_iplus1_pl )*dt
            a_i_planet = a_iplus1_pl
            v_sun += 0.5*( a_i_sun + a_iplus1_sun )*dt
            a_i_sun = a_iplus1_sun


        return x_planet,x_sun

    Planet0 = [system.initial_positions[0,0],system.initial_positions[1,0],\
    system.initial_velocities[0][0],system.initial_velocities[1][0]]

    Sun = [0,0,0,0]

    planet_orbit, sun_orbit = integrator(Planet0,Sun)
    sun_line_x = np.linspace(-7,7,100)
    sun_line_y = np.zeros_like(sun_line_x)
    sun_line_y_neg = np.zeros_like(sun_line_x)

    for i in range(len(sun_line_y)):
        sun_line_y[i] = 0.01126635
        sun_line_y_neg[i] = -0.01126635

    plt.plot(planet_orbit[:,0], planet_orbit[:,1], label="Planet 0")
    plt.plot(sun_orbit[:,0], sun_orbit[:,1], label="Sun")
    plt.plot(sun_line_x,sun_line_y , "r--", sun_line_x,sun_line_y_neg, "r--")

    plt.xlabel("x-position (AU)")
    plt.ylabel("y-position (AU)")
    plt.axis("equal")
    plt.legend()
    plt.savefig("Sun-planet-starlines.jpeg")

    plt.show()

def light_curve_our_system():
    sunmass = system.star_mass
    planetmass = system.masses[0]
    my = sunmass*planetmass/(sunmass+planetmass)
    masstot = sunmass+planetmass

    Gr = G_sol
    def integrator(planet,sun):

        def gravity_on_sun(x_sun,y_sun,x_planet,y_planet):

            r = np.asarray([x_sun-x_planet, y_sun-y_planet])
            grav = -(Gr*sunmass*planetmass)/(np.linalg.norm(r))**2
            ax = grav*r[0] /(np.linalg.norm(r))
            ay = grav*r[1]/(np.linalg.norm(r))

            return np.asarray([ax,ay],float)

        x0x,x0y,v0x,v0y = planet
        x0xsun,x0ysun,v0xsun,v0ysun = sun

        N = 1000000
        time = 0.02 #time for 20 rotations for our home planet
        #print(time)
        dt = time/N
        t = dt

        x_planet = np.zeros((N,2),float)
        x_sun = np.zeros((N,2),float)

        cm = np.array([0,0])

        x_planet[0,:] = [x0x,x0y]
        x_sun[0,:] = [x0xsun,x0ysun]

        v_planet = np.asarray([v0x,v0y],float)
        v_sun = np.asarray([v0xsun,v0ysun],float)
        F_sun0 = gravity_on_sun(x0xsun,x0ysun,x0x,x0y)
        a_i_sun = F_sun0/sunmass
        a_i_planet = -F_sun0/planetmass
        for i in range(N-1):
            t+=dt
            cm = np.asarray([(sunmass*x_sun[i,0] + planetmass*\
            x_planet[i,0])/(sunmass+planetmass),\
            (sunmass*x_sun[i,1] + planetmass*\
            x_planet[i,1])/(sunmass+planetmass)])

            x_planet[i+1,:] = x_planet[i,:]+(v_planet*dt + 0.5*a_i_planet*dt**2)-cm
            x_sun[i+1,:] = x_sun[i,:]+(v_sun*dt + 0.5*a_i_sun*dt**2)-cm
            F_sun = gravity_on_sun(x_sun[i+1,0],x_sun[i+1,1],x_planet[i+1,0],\
            x_planet[i+1,1])
            a_iplus1_sun = F_sun/sunmass
            a_iplus1_pl = -F_sun/planetmass
            v_planet += 0.5*( a_i_planet + a_iplus1_pl)*dt
            a_i_planet = a_iplus1_pl
            v_sun += 0.5*( a_i_sun + a_iplus1_sun )*dt
            a_i_sun = a_iplus1_sun

        return x_planet,x_sun , dt, N

    Planet0 = [system.initial_positions[0,0],system.initial_positions[1,0],\
    system.initial_velocities[0][0],system.initial_velocities[1][0]]

    Sun = [0,0,0,0]

    planet_orbit, sun_orbit, dt, N = integrator(Planet0,Sun)


    #From mid to
    Upper =  system.radii[0]*1000/AU + system.star_radius*1000/AU
    x_planet = np.copy(planet_orbit)
    dt = dt*yr/3600

    #Must be 1 planet radii above x axis
    condition1 = np.logical_and(np.greater_equal(x_planet[:,1],0),np.greater_equal(x_planet[:,0],0))
    condition2 = np.logical_and(np.less_equal(x_planet[:,1],Upper),np.greater_equal(x_planet[:,0],0))

    dt_in_intervall1 = np.logical_and(condition1,condition2).astype(np.float32)
    time_in_int1 = sum(dt_in_intervall1)*dt


    Upper = system.star_radius*1000/AU + system.radii[0]*1000/AU
    Under = system.star_radius*1000/AU - system.radii[0]*1000/AU

    condition3 = np.logical_and(np.greater_equal(x_planet[:,1],Under),np.greater_equal(x_planet[:,0],0))
    condition4 = np.logical_and(np.less_equal(x_planet[:,1],Upper),np.greater_equal(x_planet[:,0],0))

    dt_in_intervall2 = np.logical_and(condition3,condition4).astype(np.float32)
    time_in_int2 = sum(dt_in_intervall2*dt)

    dt_in_intervall3 = dt_in_intervall1-dt_in_intervall2

    Areal_sun = 2*np.pi*(system.star_radius*1000/AU)
    Areal_planet = 2*np.pi*(system.radii[0]*1000/AU)

    Flux_drop = (Areal_planet/Areal_sun)
    v = (Areal_planet/Areal_sun)/time_in_int2

    time = np.linspace(0,N*dt,N)
    curve = np.zeros(N)

    curve = np.where(dt_in_intervall1 == 1, dt_in_intervall1-Flux_drop, dt_in_intervall1)
    for i in range(len(dt_in_intervall2)):
        if dt_in_intervall2[i] == 1:
            curve[i] = curve[i-1] + v*dt
    curve =  np.where(curve==0,curve+1,curve)

    time = np.linspace(0,N*dt,N)
    sigma = (Areal_planet/Areal_sun)/100
    noise = np.random.normal(0,sigma, size = int(N))
    curve +=noise

    plt.title("Light-curve")
    plt.plot(time[0:int(25/dt)],curve[0:int(25/dt)])
    plt.xlabel("Time in hours")
    plt.ylabel("Flux")
    plt.show()
"""
