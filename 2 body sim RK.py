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
#from math import sin, cos, atan2
from scipy.integrate import ode

# global variables for testing, so I can see the positions results of the sim
L = np.zeros([102,4])
S = np.zeros([102,4])
    
class Body(t.Turtle):
    '''
    Defines a planetary body and its attributes. Could be planet, moon, etc.
    Inherits from Turtle for drawing purposes.
    
    Name: string which gives the body a human-readable name, such as 'Sun'
    Mass: number which gives mass of object in kg
    Pen: string which gives color of pen for animation
    Location: Array containing position in (x,y) coordinate system, in m
    Velocity: Array containing velocity in (x,y) coordinate system in m/s
    COMflag: Boolean value. True if this body is the Center Of Mass of the 
             system (i.e. the largest body)
    '''

    def __init__(self, name, mass, color, pos, vel, COM):
        t.Turtle.__init__(self)
        self.name = name
        self.mass = mass
        self.pen(pencolor=color)
        self.location = pos                         # numpy array            
        self.velocity = vel                         
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
        
        G = 6.67384 * 10**(-11)

        if self is other:
            raise ValueError("Can't attract self. This error should have been" 
                                "handled elsewhere...")

        dr = self.location - other.location         # distance between objects
        r = np.linalg.norm(dr) 
        GMr3 = -(G * other.mass) / (r**3)
        
        return GMr3


# This should maybe be a decorator function?               
def diffeqs(t, y, GMr3):
    '''
    Accepts as input y, a 4-vector of initial conditions rx, ry, vx, vy, and 
    GMr3, the 1/s^2 part that determines dv/dt. t is unused here but must be 
    passed in according to formatting needs from scipy.integrate.ode.
    '''
    return np.array([y[2], y[3], GMr3 * y[0], GMr3 * y[1]])


def solve_system(bodies, SCALE, AU):
    '''
    Solves the second order differential equation d^2r/dt^2 = -(GM/r^3)r. 
    Accepts as input the SCALE and AU variables, as they are needed here.
    '''
    
    dt = 24*3600                # Step = 1 day; in seconds
    t0 = 0
    tf = 24*3600*1001            # Total time of simulation = 100 days
    steps = tf / dt +1          # deals with indexing stupidity
    
    all_body_data = {}          # Format: [Body]: [array], steps rows, 4 col
     
    # Get the GMr3 between a body and another. Currently only works for 2 bodies
    for body in bodies:
        for other in bodies:
            if body is other:               # ignore self<>self interaction
                continue
            GMr3 = body.grav(other)
            
            # Actual diffeq evaluation
            Y = ode(diffeqs).set_integrator('dopri5')

            # Initial values. r = position, v = velocity, y0 = all, formatted
            r = body.location
            v = body.velocity
            y0 = [r[0], r[1], v[0], v[1]]
            Y.set_initial_value(y0, t0).set_f_params(GMr3)
            result_array = np.zeros([steps,4])
            
            # Loop over the time range to get solutions for all times
            while Y.successful() and Y.t <= tf:                      
                Y.integrate(Y.t + dt)
                new_row = Y.y
                i = Y.t/dt - 1
                result_array[i] = new_row
            
            time_array = Y.t
            result_array = result_array * SCALE       # don't forget to scale!
            all_body_data[body] = result_array
            
    return time_array, all_body_data


def animate(times, pos_dict):
    
    global L, S                             # for testing
    
    keys = [ key for key in pos_dict.keys() ]
    
    # Assign new variables based on which key in the dictionary has COM flag
    lg_body = keys[0] if keys[0].COMflag == True else keys[1]
    sm_body = keys[0] if keys[0].COMflag == False else keys[1]
    
    # retrieves position arrays from dictionary entries for ease of access
    lg_positions = pos_dict[lg_body]
    L = lg_positions
    sm_positions = pos_dict[sm_body]
    S = sm_positions
    
    # give planet a circle shaped turtle for visual ease
    sm_body.shape('circle')
    lg_body.penup()
    sm_body.penup()                      # for 1st step; avoid line from ctr to start
    step = 0
    
    for lg_r, sm_r in zip(lg_positions, sm_positions):
        lg_body.hideturtle()
        next_position = [lg_r[0], lg_r[1]]      # Format: rx, ry, vx, vy
        lg_body.goto(next_position)
        lg_body.dot(50)
        
        
        # Move small body to starting position without pen being down to avoid 
        # drawing a line along orbit radius
        next_position = [sm_r[0], sm_r[1]] 
        sm_body.goto(next_position)
        
        # Sets pen for small body down after body is no longer at canvas center
        if step == 2:
            sm_body.pendown()
        step += 1
    
    
def main():
    '''
    Main loop which sets the parameters and calls the simulation
    '''  
    
    # Set the mood
    t.bgcolor('black')
    
    # Assumed scale: 100 pixels = 1AU.
    AU = 149.6e9                # 149.6 billion meters is one astronomical unit
    SCALE = 250 / AU                            # 2e-9;  100 px
    
    # Declare initial positions and velocities for two bodies
    sun_loc = np.array([0,0])
    sun_vel = np.array([0,0])
    
    earth_loc = np.array([1*AU, 0])
    earth_vel = np.array([0, 29800])              # shows elliptical
    
    #mars_loc = np.array([(-227.9e9/146.6e9)*AU,0])
    #mars_vel = np.array([0,24070])
    
    # Create body objects
    sun = Body('Sun', 2e30, 'yellow', sun_loc, sun_vel, True)
    earth = Body('Earth', 5.9742e24, 'blue', earth_loc, earth_vel, False)
    #mars = Body('Mars', 6.39e23, 'red', mars_loc, mars_vel)
    
    # Do the actual work
    time, pos_data = solve_system([sun, earth], SCALE, AU)
    animate(time, pos_data)
    
if __name__ == '__main__':
    main()