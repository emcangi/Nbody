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
def diffeqs(t, y, dvx_dt, dvy_dt):
    '''
    Returns the next state in the differential equation governing motion.
    Elements are: v_x, v_y, GM/r^3 * Δx, GM/r^3 * Δy
    Vectorized; [:, 2] gets the 3rd element in each sub-array of y (r_x)
    [:, 3] gets the 4th element in each sub-array of y (r_y)
    np.sum & axis = 1 sums each sub-array and returns an array of results
    '''

    return [y[:, 2], y[:, 3], np.sum(dvx_dt, axis=1), np.sum(dvy_dt, axis=1)]


def vectorize_ic(bodies, N):
    '''
    Given the list of bodies, this function sends their initial positions and
    velocities into a 4-vector that can be used by the diffeq solver.
    '''
    i = 0
    y0 = []                         # initial conditions array

    for body in bodies:
        # Get initial values for each body, fill initial values array
        r = body.pos
        v = body.vel
        y0.append(np.array([r[0], r[1], v[0], v[1]]))
        i += 1

    return y0


def vectorize_dvdt_terms(bodies, N):
    '''
    Given the list of bodies and its length (N), collects the dv_x/dt and
    dv_y/dt terms between each body and all the others in an array. Then stores
    those arrays in two larger arrays that will serve to hold them for the
    differential equation.
    '''

    dvx_dt_all = np.ndarray(shape=(N, N-1, 1))
    dvy_dt_all = np.ndarray(shape=(N, N-1, 1))

    i = 0
    for body in bodies:
        j = 0                           # counter for filling dvx_dt, dvy_dt
        dvx_dt = np.zeros([N-1, 1])
        dvy_dt = np.zeros([N-1, 1])

        for other in bodies:
            if body is other:
                continue
            GMr3 = body.grav(other)
            dr_x = body.dist(other, 'x')
            dr_y = body.dist(other, 'y')
            dvx_dt[j] = GMr3 * dr_x
            dvy_dt[j] = GMr3 * dr_y
            j += 1
        dvx_dt_all[i] = dvx_dt
        dvy_dt_all[i] = dvy_dt
        i += 1

    return dvx_dt_all, dvy_dt_all


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
    N = len(bodies)                     # number of bodies in system

    # Make tidy array of initial conditions of all bodies
    y0 = vectorize_ic(bodies, N)

    # Make tidy arrays of the dv/dt terms
    dvx_dt, dvy_dt = vectorize_dvdt_terms(bodies, N)

    # Establish solver and integration method (Runge-Kutta, order 4)
    Y = ode(diffeqs).set_integrator('dopri5')
    Y.set_initial_value(y0, t0)
    Y.set_f_params(dvx_dt, dvy_dt)

    # Shape: Rows = time steps, columns = bodies, cells = [rx, ry, vx, vy]
    state_array = np.ndarray(shape=(t_steps, N, 4))

    # Loop over the time range to get solutions for all times
    while Y.successful() and Y.t <= tf:
        Y.integrate(Y.t + dt)
        step = Y.t/dt - 2                    # Get current step number

        # Result fmt: [[r_x1, r_y1, v_x1, v_y1], ...[r_xN, r_yN, v_xN, v_yN]].
        # Valid for this time step only!
        state_array[step] = Y.y

        # Update the body's location and velocity with new values
        for body, body_vals in zip(bodies, state_array[step]):
            body.pos = np.array([body_vals[0], body_vals[1]])
            body.vel = np.array([body_vals[1], body_vals[2]])

        # Recalculate dv/dt terms and update parameters
        dvx_dt, dvy_dt = vectorize_dvdt_terms(bodies, N)
        Y.set_f_params(dvx_dt, dvy_dt)

    time_array = Y.t

    # Convert results to AU
    state_array = state_array * SCALE

    return time_array, state_array


# ANIMATION -------------------------------------------------------------------

def run_simulation(state_array, t_steps):
    '''
    Performs animation of the simulation. Accepts as input:
    data_dict: dictionary of positions and velocities for each body over time
    t_steps: number of time steps, which will translate to number of frames

    Made with guidance from
    https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/
    '''
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    # Find the center of mass, if there is one. FIX THIS ------------------
     # FILL ME IN. IS THIS EVEN NEEDED?
    # END FIX -----------------

    # set up figure for animation. axes limits are in AU.
    fig = plt.figure()
    ax = plt.axes(xlim=(-2, 2), ylim=(-2, 2), aspect='equal')

    # Plot objects to animate--will contain vectors so we plot the ensemble
    objects, = ax.plot([], [], 'bo', ms=6)

    def init_anim():
        objects.set_data([], [])
        return objects,

    def animate(i):
        '''
        Given a frame number i (really, a given time step in our range), this
        function identifies the data associated with this time step and then
        collects all the x and y coordinates of each body into new arrays for
        easy plotting
        '''
        # Get the current row, which is one instant in time, and num of bodies
        current_state = state_array[i]
        N = len(current_state)

        # Arrays to hold data for plotting
        xdata = np.zeros([N, 1])
        ydata = np.zeros([N, 1])
        j = 0

        # collect all bodies' position data
        for bodydata in current_state:
            xdata[j] = bodydata[0]
            ydata[j] = bodydata[1]
            j += 1

        objects.set_data(xdata, ydata)

        return objects,

    anim = animation.FuncAnimation(fig, animate, frames=t_steps,
                              interval=10, blit=True, init_func=init_anim)

    anim.save('Nbody.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

    plt.show()


# MAIN ROUTINE ----------------------------------------------------------------
def main():
    '''
    Sets the parameters and calls the functions which solve the system and then
    animate it. Planet and sun sizes are not to scale.
    '''

    from math import sqrt

    planet1 = 2
    planet2 = 3

    # Initial parameters. a = semimajor axis in AU; e = eccentricity;
    # mass = mass in kg. given in arrays where first item is Mercury and last
    # item is Neptune.
    names = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn',
             'Uranus', 'Neptune']
    colors = ['r', 'y', 'b', 'r', 'r', 'm', 'c', 'b']  # matlab formatting
    size = [5, 7, 8, 5, 15, 13, 13, 12]
    axis_set = np.array([0.387098, 0.723331, 1.00000011, 1.523662,
                         5.203363, 9.537070, 19.191263, 30.068963]) * AU
    e_set = [0.20563069, 0.00677323, 0.01671022, 0.09341233,
             0.04839266, 0.0541506, 0.04716771, 0.00858587]
    mass_set = np.array([0.330, 4.87, 5.97, 0.642, 1898, 568, 86.8,
                         102]) * 10**(24)

    # Object 1
    a1 = axis_set[planet1]                      # semi-major axis, AU
    e1 = e_set[planet1]                         # eccentricity
    m1 = mass_set[planet1]
    x0_1 = a1 * (1 + e1)                        # Planet initial position
    v_y0_1 = 3e4                                # Planet initial velocity
    color1 = colors[planet1]
    name1 = names[planet1]
    size1 = size[planet1]

    # Object 2
    a2 = axis_set[planet2]                      # semi-major axis, AU
    e2 = e_set[planet2]                         # eccentricity
    m2 = mass_set[planet2]
    x0_2 = a2 * (1 + e2)                        # Planet initial position
    v_y0_2 = 3e4                                # Planet initial velocity
    color2 = colors[planet2]
    name2 = names[planet2]
    size2 = size[planet2]

    # Object 3
    m3 = 2e30
    x0_3 = 0
    v_y0_3 = 0
    color3 = 'y'
    name3 = 'Sun'
    size3 = 15

    # Declare initial positions and velocities and dt and tf
    loc1 = np.array([x0_1, 0])
    vel1 = np.array([0, v_y0_1])
    loc2 = np.array([x0_2, 0])
    vel2 = np.array([0, v_y0_2])
    loc3 = np.array([x0_3, 0])
    vel3 = np.array([0, v_y0_3])
    dt = 24*3600
    t_steps = 300

    # Create body objects
    body1 = Body(name1, m1, color1, loc1, vel1, size1, False)
    body2 = Body(name2, m2, color2, loc2, vel2, size2, False)
    body3 = Body(name3, m3, color3, loc3, vel3, size3, True)

    # Do the actual work
    time, data = solve_system([body1, body2, body3], dt, t_steps)
    run_simulation(data, t_steps)

if __name__ == '__main__':
    main()