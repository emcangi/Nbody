# -*- coding: utf-8 -*-
"""
Author: Eryn Cangi
For CIERA REU summer 2015

N-body simulation, using Newtonian mechanics and the Runge-Kutta method to
simulate interactions of N bodies. Uses matplotlib for animation.
"""

import numpy as np

G = 6.67384 * 10**(-11)
AU = 149.6e9                               # 149.6 billion meters = 1 AU

# BODY CLASS ------------------------------------------------------------------
class Body(object):
    '''
    Defines a planetary body and its attributes. Could be planet, moon, etc.

    Name: String which gives the body a human-readable name, such as 'Sun'
    Mass: Number which gives mass of object in kg
    Color: String specifying color to draw the object in, for visualization
    Location: Array containing position in (x,y) coordinate system, in m
    Velocity: Array containing velocity in (x,y) coordinate system in m/s
    Size: Integer specifying a size to draw the bodies in the animation
    COMflag: Boolean value. True if this body is the Center Of Mass of the
             system (i.e. the largest body)
    '''

    def __init__(self, name, mass, color, pos, vel, size, COM):
        self.name = name
        self.mass = mass
        self.color = color
        self.pos = pos
        self.vel = vel
        self.size = size
        self.COMflag = COM

    def __repr__(self):
        return '{name} at {pos} with speed {vel}, mass {mass}'\
               .format(name=self.name, pos=self.pos, vel=self.vel,
                       mass=self.mass)

    def grav(self, other):
        '''
        other: another Body object
        Returns the value -GM/r^3 which is needed for the diffeqs.
        '''

        if self is other:
            raise ValueError("Can't attract self. This error should have been"
                             "handled elsewhere...")

        dr = self.pos - other.pos         # distance between objects
        r = np.linalg.norm(dr)
        GMr3 = -(G * other.mass) / (r**3)

        return GMr3

    def dist(self, other, coord):
        '''
        Return signed distance between self and other along any coordinate axis
        '''
        if self is other:
            raise ValueError("Can't attract self. This error should have been"
                             "handled elsewhere...")
        if coord.lower() == 'x':
            return self.pos[0] - other.pos[0]
        elif coord.lower() == 'y':
            return self.pos[1] - other.pos[1]

# MATH ------------------------------------------------------------------------
def diffeqs(t, y, p, q):
    '''
    # FIX ME
    '''
    return # FIX ME


def solve_system(bodies, timestep, t_steps):
    '''
    Solves the second order differential equation d^2r/dt^2 = -(GM/r^3)r.
    Accepts as input a list of body objects.
    '''

    from scipy.integrate import ode

    SCALE = 1 / AU

    # Parameters
    dt = timestep
    t0 = 0
    tf = timestep * t_steps             # Total time of simulation
    #steps = tf / dt + 1
    N = len(bodies)                     # number of bodies in system

    y0 = np.zeros(N, 1)                 # initial conditions array
    p = np.zeros(N, 1)                  # array to hold each body's last2x
    q = np.zeros(N, 1)                  # array to hold each body's last2y
    all_body_data = {}                  # Format: [Body]: [steps, r, v results]

# MAKE THIS PART A SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    i = 0                               # counter for
    # For each body, solve system with reference to the other body
    for body in bodies:
        # holds (GMr3 * x dist) and (GMr3 * y dist) between body & each other
        last2x = np.zeros(N-1, 1)
        last2y = np.zeros(N-1, 1)

        j = 0                           # counter for filling previous arrays
        for other in bodies:
            if body is other:
                continue
            # INITIAL GM/r^3 and x separation, y separation btwn body, other
            GMr3 = body.grav(other)
            drx = body.dist(other, 'x')
            dry = body.dist(other, 'y')
            last2x[i] = GMr3 * drx
            last2y[i] = GMr3 * dry
            j += 1

        # when last2x, last2y complete, put in big array p or q
        p[i] = last2x
        q[i] = last2y
# END SUBROUTINE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # Get initial values for each body, fill initial values array
        r = body.pos
        v = body.vel
        y0[i] = [r[0], r[1], v[0], v[1]]
        i += 1

    # Establish solver and integration method (Runge-Kutta, order 4)
    Y = ode(diffeqs).set_integrator('dopri5')
    Y.set_initial_value(y0, t0)
    Y.set_f_params(p, q)
    result_array = np.zeros([t_steps, N])       # need N cols; 1 for each body

    # Loop over the time range to get solutions for all times
    while Y.successful() and Y.t <= tf:
        Y.integrate(Y.t + dt)
        curstep = Y.t/dt - 2                    # Get current step number
        result_array[curstep] = Y.y             # set "curstep" row to result for this body

        # Update the body's location and velocity with new values - FIX ME THIS NEEDS TO BE IN A NEW LOOP AND YOU NEED TO FIGURE OUT WHAT THE RESULTS ARE!
        body.pos = np.array([Y.y[0], Y.y[1]])
        body.vel = np.array([Y.y[2], Y.y[3]])

        # Recalculate GM/r^3 and reapply to solver - FIX ME NEED TO REDO THE LOOPS ALL OVER AGAIN OH GOD
        GMr3 = body.grav(other)
        Y.set_f_params(p, q)

    time_array = Y.t
    all_body_data[body] = result_array * SCALE

    return time_array, all_body_data


# ANIMATION -------------------------------------------------------------------

def run_simulation(time, pos_dict, t_steps):
    '''
    '''
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    # Assign variables based on which key in the dictionary has COM flag
    keys = [key for key in pos_dict.keys()]
    lg_body = keys[0] if keys[0].COMflag is True else keys[1]
    sm_body = keys[0] if keys[0].COMflag is False else keys[1]

    # retrieves position arrays from dictionary entries for ease of access
    lg_r = pos_dict[lg_body]
    sm_r = pos_dict[sm_body]

    # set up figure for animation. axes limits are in AU.
    fig = plt.figure()
    ax = plt.axes(xlim=(-2, 2), ylim=(-2, 2), aspect='equal')

    # plot element to animate (the planets)
    sun, = ax.plot([], [], lg_body.color + 'o', ms=lg_body.size)
    planet, = ax.plot([], [], sm_body.color + 'o', ms=sm_body.size)

    def init_anim():
        sun.set_data([], [])
        planet.set_data([], [])
        return sun, planet,

    def animate(i):
        sunx = lg_r[i][0]
        suny = lg_r[i][1]
        planetx = sm_r[i][0]
        planety = sm_r[i][1]
        sun.set_data(sunx, suny)
        planet.set_data(planetx, planety)
        return sun, planet,

    anim = animation.FuncAnimation(fig, animate, frames=t_steps,
                              interval=10, blit=True, init_func=init_anim)

    anim.save('2body.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    plt.show()


# MAIN ROUTINE ----------------------------------------------------------------
def main():
    '''
    Sets the parameters and calls the functions which solve the system and then
    animate it. Planet and sun sizes are not to scale.
    '''

    global time, pos_data

    from math import sqrt

    # planet selector. 0 = Mercury, 1 = Venus, ... 7 = Neptune
    planet = 2

    # Initial parameters. a = semimajor axis in AU; e = eccentricity;
    # mass = mass in kg. given in arrays where first item is Mercury and last
    # item is Neptune.
    axis_set = np.array([0.387098, 0.723331, 1.00000011, 1.523662, 5.203363,
                         9.537070, 19.191263, 30.068963]) * AU
    e_set = np.array([0.206, 0.007, 0.017, 0.093, 0.048, 0.056, 0.046, 0.010])
    mass_set = np.array([0.330, 4.87, 5.97, 0.642, 1898, 568, 86.8,
                         102]) * 10**(24)
    names = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn',
             'Uranus', 'Neptune']
    colors = ['r', 'y', 'b', 'r', 'r', 'm', 'c', 'b']  # matlab formatting
    size = [5, 7, 8, 5, 15, 13, 13, 12]

    a = axis_set[planet]
    e = e_set[planet]
    m2 = mass_set[planet]
    m1 = 2e30                                   # Larger body: the sun
    x0 = a * (1 + e)                            # initial position
    v_y0 = sqrt(G * (m1 + m2) * (2/x0 - 1/a))   # initial velocity of sm. body
    color = colors[planet]
    name = names[planet]
    size = size[planet]

    # Declare initial positions and velocities and dt and tf. Assume more
    # massive body is stationary
    sun_loc = np.array([0, 0])
    sun_vel = np.array([0, 0])
    planet_loc = np.array([x0, 0])
    planet_vel = np.array([0, v_y0])
    dt = 24*3600
    t_steps = 500

    # Create body objects
    sun = Body('Sun', m1, 'y', sun_loc, sun_vel, 20, True)
    planet = Body(name, m2, color, planet_loc, planet_vel, size, False)

    # Do the actual work
    time, pos_data = solve_system([sun, planet], dt, t_steps)
    run_simulation(time, pos_data, steps)

if __name__ == '__main__':
    main()