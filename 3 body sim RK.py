# -*- coding: utf-8 -*-
"""
3 body simulator

Author: Eryn Cangi
For CIERA REU summer 2015

3-body simulation, using Newtonian mechanics and the Runge-Kutta method to
solve the two-body problem.

Equation adapted from http://people.sc.fsu.edu/~jburkardt/m_src/three_body_simulation/three_body_simulation.html
"""

import turtle as t
import numpy as np

G = 6.67384 * 10**(-11)
AU = 149.6e9                               # 149.6 billion meters = 1 AU
SCALE = 250 / AU                           # Decrease integer for outer planets


class Body(t.Turtle):
    '''
    Defines a planetary body and its attributes. Could be planet, moon, etc.
    Inherits from Turtle for drawing purposes.

    Name: string which gives the body a human-readable name, such as 'Sun'
    Mass: number which gives mass of object in kg
    Pen: string which gives color of pen for animation
    Location: Array containing position in (x,y) coordinate system, in m
    Velocity: Array containing velocity in (x,y) coordinate system in m/s
    Size: Sets size of the planet turtle pointer. Because the adjustment is
            done with fractional values, must be divided by 13, where 13 is
            ~Earth radius. Absolutely not to scale.
    COMflag: Boolean value. True if this body is the Center Of Mass of the
             system (i.e. the largest body)
    '''

    def __init__(self, name, mass, color, pos, vel, size, COM):
        t.Turtle.__init__(self)
        self.name = name
        self.mass = mass
        self.pen(pencolor=color)
        self.location = pos             # Position is a method, so use location
        self.velocity = vel
        self.size = size
        self.COMflag = COM

    def __repr__(self):
        return '{name} at {pos} with speed {vel}, mass {mass}'\
               .format(name=self.name, pos=self.location, vel=self.velocity,
                       mass=self.mass)

    def grav(self, other):
        '''
        other: another Body object
        Returns the value -GM/r^3 which is needed for the diffeqs.
        '''

        if self is other:
            raise ValueError("Can't attract self. This error should have been"
                             "handled elsewhere...")

        dr = self.location - other.location         # distance between objects
        r = np.linalg.norm(dr)
        GMr3 = -(G * other.mass) / (r**3)

        return GMr3


def diffeqs(t, y, GMr3_0, GMr3_1, GMr3_2):
    '''
    Accepts as input y, a 4-vector of initial conditions rx, ry, vx, vy, and
    GMr3, the 1/s^2 part that determines dv/dt. t is unused here but must be
    passed in according to formatting needs from scipy.integrate.ode.
    '''
    return [y[3], y[4], GMr3_1 * (y[1] - y[5]) + GMr3_2 * (y[1] - y[9]),
            GMr3_1 * (y[2] - y[6]) + GMr3_2 * (y[2] - y[10]),
            y[7], y[8], GMr3_0 * (y[5] - y[1]) + GMr3_2 * (y[5] - y[9]),
            GMr3_0 * (y[6] - y[2]) + GMr3_2 * (y[6] - y[10]),
            y[11], y[12], GMr3_0 * (y[9] - y[1]) + GMr3_1 * (y[9] - y[5]),
            GMr3_0 * (y[10] - y[2]) + GMr3_1 * (y[10] - y[6])]


def solve_system(bodies, timestep):
    '''
    Solves the second order differential equation d^2r/dt^2 = -(GM/r^3)r.
    Accepts as input a list of body objects.
    '''

    from scipy.integrate import ode

    # Parameters
    dt = timestep                     # Step = 4 hours; in seconds
    t0 = 0
    tf = timestep*10001               # Total time of simulation = 100 days
    steps = tf / dt + 1
    ro = []
    vo = []
    GMr3 = []

    all_body_data = {}                # Format: [Body]: [steps, r, v results]

    # For each body, solve system with reference to the other body
    for body in bodies:
        for other in bodies:
            if body is other:
                continue

            # Get the INITIAL "gravitational field" (?) value from each other
            GMr3.append(body.grav(other))

            # Initial values. r = position, v = velocity; lists for other bodies
            ro.append(other.location)
            vo.append(other.location)

        # I think all this stuff has to happen outside of the second loop actually
        # Establish solver and integration method (Runge-Kutta, order 4)
        Y = ode(diffeqs).set_integrator('dopri5')
        r = body.location
        v = body.velocity
        y0 = [r[0], r[1], v[0], v[1], ro[0], ro[1], vo[0], vo[1], ro[2],
              ro[3], vo[2], vo[3]]
        Y.set_initial_value(y0, t0)
        Y.set_f_params(GMr3)
        result_array = np.zeros([steps, 12])

        # Loop over the time range to get solutions for all times
        while Y.successful() and Y.t <= tf:
            Y.integrate(Y.t + dt)
            i = Y.t/dt - 1                        # Get current step number
            result_array[i] = Y.y

            # Update the body's location and velocity with new values
            body.location = np.array([Y.y[0], Y.y[1]])
            body.velocity = np.array([Y.y[2], Y.y[3]])

            # HOW TO UPDATE THE POSITIONS, VELOCITIES OF THE OTHER BODIES?
            # ANOTHER LOOP?

            # Recalculate GM/r^3 and reapply to solver
            GMr3 = body.grav(other)
            Y.set_f_params(GMr3, r_other)

            time_array = Y.t
            all_body_data[body] = result_array * SCALE

    return time_array, all_body_data


def animate(times, pos_dict):
    '''
    Takes data generated by the system solver and animates it using Turtle
    graphics.
    '''

    # Assign variables based on which key in the dictionary has COM flag
    keys = [key for key in pos_dict.keys()]
    lg_body = keys[0] if keys[0].COMflag is True else keys[1]
    sm_body = keys[0] if keys[0].COMflag is False else keys[1]

    # retrieves position arrays from dictionary entries for ease of access
    lg_positions = pos_dict[lg_body]
    sm_positions = pos_dict[sm_body]

    # give planet a circle shaped turtle for familiar visuals
    sm_body.shape('circle')
    sm_body.resizemode("user")
    n = sm_body.size
    sm_body.shapesize(n, n, n)
    lg_body.penup()
    sm_body.penup()

    for lg_r, sm_r in zip(lg_positions, sm_positions):
        lg_body.hideturtle()
        next_position = [lg_r[0], lg_r[1]]      # Format: rx, ry, vx, vy
        lg_body.goto(next_position)
        lg_body.dot(lg_body.size)

        # Move small body to orbit path w/ pen up to avoid line on orbit radius
        next_position = [sm_r[0], sm_r[1]]
        sm_body.goto(next_position)
        sm_body.pendown()


def main():
    '''
    Sets the parameters and calls the functions which solve the system and then
    animate it. Planet and sun sizes are not to scale.
    '''

    from math import sqrt

    # Set the mood
    t.bgcolor('black')

    # planet selector. 0 = Mercury, 1 = Venus, ... 7 = Neptune
    planet = 7

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
    colors = ['red', 'yellow', 'blue', 'red', 'red', 'orange', 'cyan', 'blue']
    pensize = [1/13, 3/13, 4/13, 26/13, 2, 26/13, 20/13, 19/13]

    a = axis_set[planet]
    e = e_set[planet]
    m2 = mass_set[planet]
    m1 = 2e30                                   # Larger body: the sun
    x0 = a * (1 + e)                            # initial position
    v_y0 = sqrt(G * (m1 + m2) * (2/x0 - 1/a))   # initial velocity of sm. body
    color = colors[planet]
    name = names[planet]
    size = pensize[planet]

    # Declare initial positions and velocities for two bodies. Assume more
    # massive body is stationary
    sun_loc = np.array([0, 0])
    sun_vel = np.array([0, 0])
    planet_loc = np.array([x0, 0])
    planet_vel = np.array([0, v_y0])

    # Create body objects
    sun = Body('Sun', m1, 'yellow', sun_loc, sun_vel, 40, True)
    planet = Body(name, m2, color, planet_loc, planet_vel, size, False)

    # Do the actual work
    time, pos_data = solve_system([sun, planet], 24*3600)
    animate(time, pos_data)

if __name__ == '__main__':
    main()