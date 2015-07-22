# -*- coding: utf-8 -*-
"""
2 body simulator

Author: Eryn Cangi
For CIERA REU summer 2015

2-body simulation, using Newtonian mechanics and the Runge-Kutta method to 
solve the two-body problem. 

Heavily adapted from http://fiftyexamples.readthedocs.org/en/latest/gravity.html
License is MIT license, free to use and modify as desired for any use.
"""

import turtle as t
import numpy as np

G = 6.67384 * 10**(-11)
AU = 149.6e9                                # 149.6 billion meters = 1 AU
SCALE = 250 / AU
    
class Body(t.Turtle):
    '''
    Defines a planetary body and its attributes. Could be planet, moon, etc.
    Inherits from Turtle for drawing purposes.
    
    Name: string which gives the body a human-readable name, such as 'Sun'
    Mass: number which gives mass of object in kg
    Pen: string which gives color of pen for animation
    Location: Array containing position in (x,y) coordinate system, in m
    Velocity: Array containing velocity in (x,y) coordinate system in m/s
    Size: Integer value roughly corresponding to diameter in km, divided by 13,
            where 13 is ~Earth radius. Division is required because pen size 
            for small body can only be set as percentage of normal size.
    COMflag: Boolean value. True if this body is the Center Of Mass of the 
             system (i.e. the largest body)
    '''

    def __init__(self, name, mass, color, pos, vel, size, COM):
        t.Turtle.__init__(self)
        self.name = name
        self.mass = mass
        self.pen(pencolor=color)
        self.location = pos                         # numpy array            
        self.velocity = vel
        self.size = size                    
        self.COMflag = COM

    def __repr__(self):
        return '{name} at {pos} with speed {vel}, mass {mass}'\
        .format(name = self.name, pos = self.location, vel = self.velocity, 
                mass = self.mass)
        
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
        GMr3 = -(G * (other.mass + self.mass)) / (r**3)
        
        return GMr3


# This should maybe be a decorator function?               
def diffeqs(t, y, GMr3):
    '''
    Accepts as input y, a 4-vector of initial conditions rx, ry, vx, vy, and 
    GMr3, the 1/s^2 part that determines dv/dt. t is unused here but must be 
    passed in according to formatting needs from scipy.integrate.ode.
    '''
    return np.array([y[2], y[3], GMr3 * y[0], GMr3 * y[1]])


def solve_system(bodies):
    '''
    Solves the second order differential equation d^2r/dt^2 = -(GM/r^3)r. 
    Accepts as input a list of body objects.
    '''

    from scipy.integrate import ode

    # Parameters    
    dt = 24*360                    # Step = 1 day; in seconds
    t0 = 0
    tf = 24*360*1001               # Total time of simulation = 100 days
    steps = tf / dt +1              # deals with indexing stupidity
    
    all_body_data = {}              # Format: [Body]: [timesteps, r and v results]
    
    # For each body, solve system with reference to the other body
    for body in bodies:
        for other in bodies:
            if body is other:
                continue
            
            # Get the INITIAL "gravitational field" (?) value
            GMr3 = body.grav(other) 
                
            # Establish solver and integration method (Runge-Kutta, order 4)
            Y = ode(diffeqs).set_integrator('dopri5')

            # Initial values. r = position, v = velocity, y0 = formatted vector
            r = body.location
            v = body.velocity
            y0 = [r[0], r[1], v[0], v[1]]
            Y.set_initial_value(y0, t0)
            Y.set_f_params(GMr3)
            result_array = np.zeros([steps,4])
            
            # Loop over the time range to get solutions for all times
            while Y.successful() and Y.t <= tf:                      
                Y.integrate(Y.t + dt)
                i = Y.t/dt - 1                        # Get current step number
                result_array[i] = Y.y   
                
                # Update the body's location and velocity with new values
                body.location = np.array([Y.y[0], Y.y[1]])
                body.velocity = np.array([Y.y[2], Y.y[3]])
                
                # Recalculate the gravitational "Field" (?) and reapply to solver
                GMr3 = body.grav(other)
                Y.set_f_params(GMr3)
                
            time_array = Y.t
            all_body_data[body] = result_array * SCALE
            
    return time_array, all_body_data


def animate(times, pos_dict):
    '''
    Takes data generated by the system solver and animates it using Turtle 
    graphics. 
    '''
    
    # Assign variables based on which key in the dictionary has COM flag
    keys = [ key for key in pos_dict.keys() ]
    lg_body = keys[0] if keys[0].COMflag == True else keys[1]
    sm_body = keys[0] if keys[0].COMflag == False else keys[1]
    
    # retrieves position arrays from dictionary entries for ease of access
    lg_positions = pos_dict[lg_body]
    sm_positions = pos_dict[sm_body]
    
    # give planet a circle shaped turtle for familiar visuals
    sm_body.shape('circle')
    sm_body.resizemode("user")
    n = sm_body.size
    sm_body.shapesize(n,n,n)
    lg_body.penup()
    sm_body.penup()
    
    for lg_r, sm_r in zip(lg_positions, sm_positions):
        lg_body.hideturtle()
        next_position = [lg_r[0], lg_r[1]]      # Format: rx, ry, vx, vy
        lg_body.goto(next_position)
        lg_body.dot(35)
        
        # Move small body to orbit path w/ pen up to avoid line on orbit radius
        next_position = [sm_r[0], sm_r[1]] 
        sm_body.goto(next_position)
        sm_body.pendown()
        
    
def main():
    '''
    Sets the parameters and calls the functions which solve the system and then
    animate it.
    '''  
    
    from math import sqrt
    
    # Set the mood
    t.bgcolor('black')
    
    # planet selector. 0 = Mercury, 1 = Venus, ... 7 = Neptune
    planet = 0
    
    # Initial parameters. a = semimajor axis in AU; e = eccentricity; 
    # mass = mass in kg. given in arrays where first item is Mercury and last 
    # item is Neptune.
    axis_set = np.array([0.387098, 0.723331, 1.00000011, 1.523662, 5.203363, 
                          9.537070, 19.191263, 30.068963])* AU                              
    e_set = np.array([0.206, 0.007, 0.017, 0.093, 0.048, 0.056, 0.046, 0.010])
    mass_set = np.array([0.330, 4.87, 5.97, 0.642, 1898, 568, 86.8, 
                         102]) * 10**(24)
    names = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 
                'Uranus', 'Neptune']
    colors = ['red', 'yellow', 'blue', 'red', 'red', 'orange', 'cyan', 'blue']
    pensize = [5/13, 12/13, 1, 7/13, 143/13, 120/13, 51/13, 50/13] 
    
    a = axis_set[planet]
    e = e_set[planet]
    m2 = mass_set[planet]
    m1 = 2e30                                   # Larger body                 
    x0 = a * (1 + e)                            # initial position
    v_y0 = sqrt(G * (m1 + m2) * (2/x0 - 1/a))   # initial velocity of sm. body
    color = colors[planet]
    name = names[planet]
    size = pensize[planet]
    
    # Declare initial positions and velocities for two bodies. Assume more 
    # massive body is stationary
    sun_loc = np.array([0,0])
    sun_vel = np.array([0,0])
    planet_loc = np.array([x0,0])
    planet_vel = np.array([0,v_y0])
    
    # Create body objects
    sun = Body('Sun', m1, 'yellow', sun_loc, sun_vel, 200, True)
    planet = Body(name, m2, color, planet_loc, planet_vel, size, False)
    
    # Do the actual work
    time, pos_data = solve_system([sun, planet])
    animate(time, pos_data)
    
if __name__ == '__main__':
    main()