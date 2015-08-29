function [vxf1, vyf1, vxf2, vyf2, fail] = collision(vxi1, vyi1, vxi2, vyi2)
% Collision handler for N-body simulation. Solution for final velocities is
% horrifyingly ugly and comes from a Maple script. x represents the
% percentage of total energy retained after the collision, and f represents
% the fraction of the energy assigned to body 1. It is restricted to the
% range 0.20 to 0.80 to avoid violating conservation of momentum.

x = 0.9;
f = (0.60 - 0.40) * rand(1) + 0.30;
fail = 0;


vxf1 = -(((2 * vyi1 * x * vxi1 ^ 2 * f) + (2 * vyi1 * x * vxi2 ^ 2 * f) + (2 * vyi1 * x * vyi2 ^ 2 * f) + (2 * vyi2 * x * vxi1 ^ 2 * f) + (2 * vyi2 * x * vxi2 ^ 2 * f) + (2 * vyi2 * x * vyi1 ^ 2 * f) + (2 * vyi2 ^ 3 * x * f) + (2 * vyi1 ^ 3 * x * f) + (2 * vyi1 * vxi1 * vxi2) - (vyi1 * x * vxi1 ^ 2) - (vyi1 * x * vxi2 ^ 2) - (vyi1 * x * vyi2 ^ 2) + (2 * vyi2 * vxi1 * vxi2) - (vyi2 * x * vxi1 ^ 2) - (vyi2 * x * vyi1 ^ 2) - (vyi2 * x * vxi2 ^ 2) + sqrt((-16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f ^ 2 + 16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f - 2 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 - 2 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 - 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f ^ 2 + 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f - 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f ^ 2 + 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi2 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi2 ^ 2 + 12 * vxi1 ^ 4 * x ^ 2 * f * vxi2 ^ 2 - 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 12 * vxi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi2 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 + 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f + 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 12 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi1 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 16 * vyi1 * vxi2 ^ 3 * vyi2 * vxi1 + 4 * vyi1 * vxi2 ^ 4 * vyi2 * x - 8 * vyi1 * vyi2 ^ 3 * vxi1 * vxi2 + 16 * vyi2 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi1 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi2 ^ 2 + 4 * vyi1 ^ 4 * x * vxi1 * vxi2 + 4 * vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 * f + 4 * vyi2 ^ 3 * x * vyi1 * vxi1 ^ 2 + 4 * vyi2 ^ 3 * x * vyi1 * vxi2 ^ 2 + 4 * vyi2 ^ 4 * x * vxi1 * vxi2 + 4 * vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 * f - 8 * vyi1 ^ 3 * vyi2 * vxi1 * vxi2 - 12 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 * vxi2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi1 ^ 2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi2 ^ 2 - 24 * vyi1 * vxi1 ^ 2 * vyi2 * vxi2 ^ 2 + 16 * vyi1 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 16 * vyi1 * vxi1 ^ 3 * vyi2 * vxi2 + 4 * vyi1 * vxi1 ^ 4 * vyi2 * x + 2 * vxi1 ^ 6 * x - vxi1 ^ 6 * x ^ 2 - 4 * vxi1 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi1 ^ 6 * x ^ 2 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 + 2 * vyi1 ^ 4 * x * vxi1 ^ 2 + 2 * vyi1 ^ 4 * x * vxi2 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 2 * vyi2 ^ 4 * x * vxi1 ^ 2 + 2 * vyi2 ^ 4 * x * vxi2 ^ 2 - vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 - 2 * vxi2 * vxi1 ^ 5 * x ^ 2 - 2 * vxi2 ^ 5 * vxi1 * x ^ 2 - 4 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 + 14 * vxi1 ^ 2 * vxi2 ^ 4 * x - 3 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 - 3 * vxi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 8 * vxi1 ^ 5 * vxi2 * x + 16 * vxi1 ^ 3 * vxi2 ^ 3 * x + 14 * vxi1 ^ 4 * x * vxi2 ^ 2 - 4 * vxi2 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi2 ^ 6 * x ^ 2 * f + 8 * vxi2 ^ 5 * vxi1 * x - 2 * vyi2 ^ 4 * vxi1 * vxi2 - 2 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi1 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi2 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi2 ^ 2 - 12 * vyi1 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 4 * vyi1 * vxi1 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi1 ^ 2 * vxi1 ^ 4 * x - 4 * vyi1 * vxi1 ^ 2 * vyi2 ^ 3 - 4 * vyi1 * vxi2 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi1 ^ 2 * vxi2 ^ 4 * x - 4 * vyi1 * vxi2 ^ 2 * vyi2 ^ 3 - 12 * vyi2 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 8 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi2 ^ 2 * vxi1 ^ 4 * x - 8 * vyi2 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi2 ^ 2 * vxi2 ^ 4 * x - 2 * vyi1 ^ 4 * vxi1 * vxi2 - vxi1 ^ 6 - vxi2 ^ 6 - 2 * vyi2 ^ 2 * vxi1 ^ 4 - vyi2 ^ 4 * vxi1 ^ 2 - 2 * vyi2 ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * vxi2 ^ 2 - 15 * vxi1 ^ 4 * vxi2 ^ 2 - 6 * vxi1 ^ 5 * vxi2 - 20 * vxi1 ^ 3 * vxi2 ^ 3 - 15 * vxi1 ^ 2 * vxi2 ^ 4 - 6 * vxi2 ^ 5 * vxi1 + 2 * vxi2 ^ 6 * x - vxi2 ^ 6 * x ^ 2 - 2 * vyi1 ^ 2 * vxi1 ^ 4 - vyi1 ^ 4 * vxi1 ^ 2 - 2 * vyi1 ^ 2 * vxi2 ^ 4 - vyi1 ^ 4 * vxi2 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi1 ^ 2 - 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f - 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f - 4 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 2 * vyi2 ^ 2 - 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vyi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * f ^ 2 * vyi2 ^ 2 + 8 * vyi1 * vxi1 ^ 3 * vxi2 * vyi2 * x + 8 * vyi1 * vxi1 * vxi2 ^ 3 * vyi2 * x + 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 * f + 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 * vxi1 ^ 2 * vyi2 * x * vxi2 ^ 2 + 8 * vyi1 ^ 2 * vxi1 * vxi2 * x * vyi2 ^ 2 + 8 * vyi1 ^ 3 * x * vyi2 * vxi1 * vxi2 + 8 * vyi2 ^ 3 * x * vyi1 * vxi1 * vxi2 - 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f * vyi1 ^ 2)) + (vyi1 ^ 3) + (vyi2 ^ 3) - (vyi1 ^ 3 * x) - (vyi2 ^ 3 * x) + (3 * vyi1 ^ 2 * vyi2) + (vyi1 * vxi1 ^ 2) + (vyi1 * vxi2 ^ 2) + (3 * vyi1 * vyi2 ^ 2) + (vyi2 * vxi1 ^ 2) + (vyi2 * vxi2 ^ 2)) / (vyi2 ^ 2 + vxi1 ^ 2 + 2 * vyi2 * vyi1 + vxi2 ^ 2 + vyi1 ^ 2 + 2 * vxi2 * vxi1) * vyi1 + ((2 * vyi1 * x * vxi1 ^ 2 * f) + (2 * vyi1 * x * vxi2 ^ 2 * f) + (2 * vyi1 * x * vyi2 ^ 2 * f) + (2 * vyi2 * x * vxi1 ^ 2 * f) + (2 * vyi2 * x * vxi2 ^ 2 * f) + (2 * vyi2 * x * vyi1 ^ 2 * f) + (2 * vyi2 ^ 3 * x * f) + (2 * vyi1 ^ 3 * x * f) + (2 * vyi1 * vxi1 * vxi2) - (vyi1 * x * vxi1 ^ 2) - (vyi1 * x * vxi2 ^ 2) - (vyi1 * x * vyi2 ^ 2) + (2 * vyi2 * vxi1 * vxi2) - (vyi2 * x * vxi1 ^ 2) - (vyi2 * x * vyi1 ^ 2) - (vyi2 * x * vxi2 ^ 2) + sqrt((-16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f ^ 2 + 16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f - 2 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 - 2 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 - 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f ^ 2 + 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f - 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f ^ 2 + 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi2 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi2 ^ 2 + 12 * vxi1 ^ 4 * x ^ 2 * f * vxi2 ^ 2 - 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 12 * vxi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi2 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 + 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f + 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 12 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi1 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 16 * vyi1 * vxi2 ^ 3 * vyi2 * vxi1 + 4 * vyi1 * vxi2 ^ 4 * vyi2 * x - 8 * vyi1 * vyi2 ^ 3 * vxi1 * vxi2 + 16 * vyi2 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi1 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi2 ^ 2 + 4 * vyi1 ^ 4 * x * vxi1 * vxi2 + 4 * vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 * f + 4 * vyi2 ^ 3 * x * vyi1 * vxi1 ^ 2 + 4 * vyi2 ^ 3 * x * vyi1 * vxi2 ^ 2 + 4 * vyi2 ^ 4 * x * vxi1 * vxi2 + 4 * vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 * f - 8 * vyi1 ^ 3 * vyi2 * vxi1 * vxi2 - 12 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 * vxi2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi1 ^ 2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi2 ^ 2 - 24 * vyi1 * vxi1 ^ 2 * vyi2 * vxi2 ^ 2 + 16 * vyi1 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 16 * vyi1 * vxi1 ^ 3 * vyi2 * vxi2 + 4 * vyi1 * vxi1 ^ 4 * vyi2 * x + 2 * vxi1 ^ 6 * x - vxi1 ^ 6 * x ^ 2 - 4 * vxi1 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi1 ^ 6 * x ^ 2 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 + 2 * vyi1 ^ 4 * x * vxi1 ^ 2 + 2 * vyi1 ^ 4 * x * vxi2 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 2 * vyi2 ^ 4 * x * vxi1 ^ 2 + 2 * vyi2 ^ 4 * x * vxi2 ^ 2 - vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 - 2 * vxi2 * vxi1 ^ 5 * x ^ 2 - 2 * vxi2 ^ 5 * vxi1 * x ^ 2 - 4 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 + 14 * vxi1 ^ 2 * vxi2 ^ 4 * x - 3 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 - 3 * vxi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 8 * vxi1 ^ 5 * vxi2 * x + 16 * vxi1 ^ 3 * vxi2 ^ 3 * x + 14 * vxi1 ^ 4 * x * vxi2 ^ 2 - 4 * vxi2 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi2 ^ 6 * x ^ 2 * f + 8 * vxi2 ^ 5 * vxi1 * x - 2 * vyi2 ^ 4 * vxi1 * vxi2 - 2 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi1 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi2 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi2 ^ 2 - 12 * vyi1 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 4 * vyi1 * vxi1 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi1 ^ 2 * vxi1 ^ 4 * x - 4 * vyi1 * vxi1 ^ 2 * vyi2 ^ 3 - 4 * vyi1 * vxi2 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi1 ^ 2 * vxi2 ^ 4 * x - 4 * vyi1 * vxi2 ^ 2 * vyi2 ^ 3 - 12 * vyi2 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 8 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi2 ^ 2 * vxi1 ^ 4 * x - 8 * vyi2 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi2 ^ 2 * vxi2 ^ 4 * x - 2 * vyi1 ^ 4 * vxi1 * vxi2 - vxi1 ^ 6 - vxi2 ^ 6 - 2 * vyi2 ^ 2 * vxi1 ^ 4 - vyi2 ^ 4 * vxi1 ^ 2 - 2 * vyi2 ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * vxi2 ^ 2 - 15 * vxi1 ^ 4 * vxi2 ^ 2 - 6 * vxi1 ^ 5 * vxi2 - 20 * vxi1 ^ 3 * vxi2 ^ 3 - 15 * vxi1 ^ 2 * vxi2 ^ 4 - 6 * vxi2 ^ 5 * vxi1 + 2 * vxi2 ^ 6 * x - vxi2 ^ 6 * x ^ 2 - 2 * vyi1 ^ 2 * vxi1 ^ 4 - vyi1 ^ 4 * vxi1 ^ 2 - 2 * vyi1 ^ 2 * vxi2 ^ 4 - vyi1 ^ 4 * vxi2 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi1 ^ 2 - 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f - 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f - 4 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 2 * vyi2 ^ 2 - 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vyi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * f ^ 2 * vyi2 ^ 2 + 8 * vyi1 * vxi1 ^ 3 * vxi2 * vyi2 * x + 8 * vyi1 * vxi1 * vxi2 ^ 3 * vyi2 * x + 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 * f + 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 * vxi1 ^ 2 * vyi2 * x * vxi2 ^ 2 + 8 * vyi1 ^ 2 * vxi1 * vxi2 * x * vyi2 ^ 2 + 8 * vyi1 ^ 3 * x * vyi2 * vxi1 * vxi2 + 8 * vyi2 ^ 3 * x * vyi1 * vxi1 * vxi2 - 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f * vyi1 ^ 2)) + (vyi1 ^ 3) + (vyi2 ^ 3) - (vyi1 ^ 3 * x) - (vyi2 ^ 3 * x) + (3 * vyi1 ^ 2 * vyi2) + (vyi1 * vxi1 ^ 2) + (vyi1 * vxi2 ^ 2) + (3 * vyi1 * vyi2 ^ 2) + (vyi2 * vxi1 ^ 2) + (vyi2 * vxi2 ^ 2)) / (vyi2 ^ 2 + vxi1 ^ 2 + 2 * vyi2 * vyi1 + vxi2 ^ 2 + vyi1 ^ 2 + 2 * vxi2 * vxi1) * vyi2 - (2 * vyi2 * vyi1) - (2 * vxi2 * vxi1) - (vxi1 ^ 2) - (vyi1 ^ 2) - (vxi2 ^ 2) - (vyi2 ^ 2) - (2 * x * vxi1 ^ 2 * f) - (2 * x * vxi2 ^ 2 * f) - (2 * x * vyi2 ^ 2 * f) - (2 * x * vyi1 ^ 2 * f) + (x * vxi1 ^ 2) + (vyi1 ^ 2 * x) + (vxi2 ^ 2 * x) + (vyi2 ^ 2 * x)) / (vxi1 + vxi2) / 0.2e1;
vyf1 = ((2 * vyi1 * x * vxi1 ^ 2 * f) + (2 * vyi1 * x * vxi2 ^ 2 * f) + (2 * vyi1 * x * vyi2 ^ 2 * f) + (2 * vyi2 * x * vxi1 ^ 2 * f) + (2 * vyi2 * x * vxi2 ^ 2 * f) + (2 * vyi2 * x * vyi1 ^ 2 * f) + (2 * vyi2 ^ 3 * x * f) + (2 * vyi1 ^ 3 * x * f) + (2 * vyi1 * vxi1 * vxi2) - (vyi1 * x * vxi1 ^ 2) - (vyi1 * x * vxi2 ^ 2) - (vyi1 * x * vyi2 ^ 2) + (2 * vyi2 * vxi1 * vxi2) - (vyi2 * x * vxi1 ^ 2) - (vyi2 * x * vyi1 ^ 2) - (vyi2 * x * vxi2 ^ 2) + sqrt((-16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f ^ 2 + 16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f - 2 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 - 2 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 - 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f ^ 2 + 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f - 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f ^ 2 + 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi2 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi2 ^ 2 + 12 * vxi1 ^ 4 * x ^ 2 * f * vxi2 ^ 2 - 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 12 * vxi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi2 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 + 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f + 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 12 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi1 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 16 * vyi1 * vxi2 ^ 3 * vyi2 * vxi1 + 4 * vyi1 * vxi2 ^ 4 * vyi2 * x - 8 * vyi1 * vyi2 ^ 3 * vxi1 * vxi2 + 16 * vyi2 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi1 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi2 ^ 2 + 4 * vyi1 ^ 4 * x * vxi1 * vxi2 + 4 * vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 * f + 4 * vyi2 ^ 3 * x * vyi1 * vxi1 ^ 2 + 4 * vyi2 ^ 3 * x * vyi1 * vxi2 ^ 2 + 4 * vyi2 ^ 4 * x * vxi1 * vxi2 + 4 * vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 * f - 8 * vyi1 ^ 3 * vyi2 * vxi1 * vxi2 - 12 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 * vxi2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi1 ^ 2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi2 ^ 2 - 24 * vyi1 * vxi1 ^ 2 * vyi2 * vxi2 ^ 2 + 16 * vyi1 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 16 * vyi1 * vxi1 ^ 3 * vyi2 * vxi2 + 4 * vyi1 * vxi1 ^ 4 * vyi2 * x + 2 * vxi1 ^ 6 * x - vxi1 ^ 6 * x ^ 2 - 4 * vxi1 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi1 ^ 6 * x ^ 2 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 + 2 * vyi1 ^ 4 * x * vxi1 ^ 2 + 2 * vyi1 ^ 4 * x * vxi2 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 2 * vyi2 ^ 4 * x * vxi1 ^ 2 + 2 * vyi2 ^ 4 * x * vxi2 ^ 2 - vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 - 2 * vxi2 * vxi1 ^ 5 * x ^ 2 - 2 * vxi2 ^ 5 * vxi1 * x ^ 2 - 4 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 + 14 * vxi1 ^ 2 * vxi2 ^ 4 * x - 3 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 - 3 * vxi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 8 * vxi1 ^ 5 * vxi2 * x + 16 * vxi1 ^ 3 * vxi2 ^ 3 * x + 14 * vxi1 ^ 4 * x * vxi2 ^ 2 - 4 * vxi2 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi2 ^ 6 * x ^ 2 * f + 8 * vxi2 ^ 5 * vxi1 * x - 2 * vyi2 ^ 4 * vxi1 * vxi2 - 2 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi1 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi2 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi2 ^ 2 - 12 * vyi1 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 4 * vyi1 * vxi1 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi1 ^ 2 * vxi1 ^ 4 * x - 4 * vyi1 * vxi1 ^ 2 * vyi2 ^ 3 - 4 * vyi1 * vxi2 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi1 ^ 2 * vxi2 ^ 4 * x - 4 * vyi1 * vxi2 ^ 2 * vyi2 ^ 3 - 12 * vyi2 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 8 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi2 ^ 2 * vxi1 ^ 4 * x - 8 * vyi2 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi2 ^ 2 * vxi2 ^ 4 * x - 2 * vyi1 ^ 4 * vxi1 * vxi2 - vxi1 ^ 6 - vxi2 ^ 6 - 2 * vyi2 ^ 2 * vxi1 ^ 4 - vyi2 ^ 4 * vxi1 ^ 2 - 2 * vyi2 ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * vxi2 ^ 2 - 15 * vxi1 ^ 4 * vxi2 ^ 2 - 6 * vxi1 ^ 5 * vxi2 - 20 * vxi1 ^ 3 * vxi2 ^ 3 - 15 * vxi1 ^ 2 * vxi2 ^ 4 - 6 * vxi2 ^ 5 * vxi1 + 2 * vxi2 ^ 6 * x - vxi2 ^ 6 * x ^ 2 - 2 * vyi1 ^ 2 * vxi1 ^ 4 - vyi1 ^ 4 * vxi1 ^ 2 - 2 * vyi1 ^ 2 * vxi2 ^ 4 - vyi1 ^ 4 * vxi2 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi1 ^ 2 - 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f - 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f - 4 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 2 * vyi2 ^ 2 - 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vyi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * f ^ 2 * vyi2 ^ 2 + 8 * vyi1 * vxi1 ^ 3 * vxi2 * vyi2 * x + 8 * vyi1 * vxi1 * vxi2 ^ 3 * vyi2 * x + 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 * f + 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 * vxi1 ^ 2 * vyi2 * x * vxi2 ^ 2 + 8 * vyi1 ^ 2 * vxi1 * vxi2 * x * vyi2 ^ 2 + 8 * vyi1 ^ 3 * x * vyi2 * vxi1 * vxi2 + 8 * vyi2 ^ 3 * x * vyi1 * vxi1 * vxi2 - 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f * vyi1 ^ 2)) + (vyi1 ^ 3) + (vyi2 ^ 3) - (vyi1 ^ 3 * x) - (vyi2 ^ 3 * x) + (3 * vyi1 ^ 2 * vyi2) + (vyi1 * vxi1 ^ 2) + (vyi1 * vxi2 ^ 2) + (3 * vyi1 * vyi2 ^ 2) + (vyi2 * vxi1 ^ 2) + (vyi2 * vxi2 ^ 2)) / (vyi2 ^ 2 + vxi1 ^ 2 + 2 * vyi2 * vyi1 + vxi2 ^ 2 + vyi1 ^ 2 + 2 * vxi2 * vxi1) / 0.2e1;

vxf2 = (((2 * vyi1 * x * vxi1 ^ 2 * f) + (2 * vyi1 * x * vxi2 ^ 2 * f) + (2 * vyi1 * x * vyi2 ^ 2 * f) + (2 * vyi2 * x * vxi1 ^ 2 * f) + (2 * vyi2 * x * vxi2 ^ 2 * f) + (2 * vyi2 * x * vyi1 ^ 2 * f) + (2 * vyi2 ^ 3 * x * f) + (2 * vyi1 ^ 3 * x * f) + (2 * vyi1 * vxi1 * vxi2) - (vyi1 * x * vxi1 ^ 2) - (vyi1 * x * vxi2 ^ 2) - (vyi1 * x * vyi2 ^ 2) + (2 * vyi2 * vxi1 * vxi2) - (vyi2 * x * vxi1 ^ 2) - (vyi2 * x * vyi1 ^ 2) - (vyi2 * x * vxi2 ^ 2) + sqrt((-16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f ^ 2 + 16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f - 2 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 - 2 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 - 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f ^ 2 + 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f - 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f ^ 2 + 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi2 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi2 ^ 2 + 12 * vxi1 ^ 4 * x ^ 2 * f * vxi2 ^ 2 - 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 12 * vxi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi2 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 + 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f + 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 12 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi1 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 16 * vyi1 * vxi2 ^ 3 * vyi2 * vxi1 + 4 * vyi1 * vxi2 ^ 4 * vyi2 * x - 8 * vyi1 * vyi2 ^ 3 * vxi1 * vxi2 + 16 * vyi2 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi1 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi2 ^ 2 + 4 * vyi1 ^ 4 * x * vxi1 * vxi2 + 4 * vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 * f + 4 * vyi2 ^ 3 * x * vyi1 * vxi1 ^ 2 + 4 * vyi2 ^ 3 * x * vyi1 * vxi2 ^ 2 + 4 * vyi2 ^ 4 * x * vxi1 * vxi2 + 4 * vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 * f - 8 * vyi1 ^ 3 * vyi2 * vxi1 * vxi2 - 12 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 * vxi2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi1 ^ 2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi2 ^ 2 - 24 * vyi1 * vxi1 ^ 2 * vyi2 * vxi2 ^ 2 + 16 * vyi1 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 16 * vyi1 * vxi1 ^ 3 * vyi2 * vxi2 + 4 * vyi1 * vxi1 ^ 4 * vyi2 * x + 2 * vxi1 ^ 6 * x - vxi1 ^ 6 * x ^ 2 - 4 * vxi1 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi1 ^ 6 * x ^ 2 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 + 2 * vyi1 ^ 4 * x * vxi1 ^ 2 + 2 * vyi1 ^ 4 * x * vxi2 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 2 * vyi2 ^ 4 * x * vxi1 ^ 2 + 2 * vyi2 ^ 4 * x * vxi2 ^ 2 - vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 - 2 * vxi2 * vxi1 ^ 5 * x ^ 2 - 2 * vxi2 ^ 5 * vxi1 * x ^ 2 - 4 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 + 14 * vxi1 ^ 2 * vxi2 ^ 4 * x - 3 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 - 3 * vxi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 8 * vxi1 ^ 5 * vxi2 * x + 16 * vxi1 ^ 3 * vxi2 ^ 3 * x + 14 * vxi1 ^ 4 * x * vxi2 ^ 2 - 4 * vxi2 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi2 ^ 6 * x ^ 2 * f + 8 * vxi2 ^ 5 * vxi1 * x - 2 * vyi2 ^ 4 * vxi1 * vxi2 - 2 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi1 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi2 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi2 ^ 2 - 12 * vyi1 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 4 * vyi1 * vxi1 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi1 ^ 2 * vxi1 ^ 4 * x - 4 * vyi1 * vxi1 ^ 2 * vyi2 ^ 3 - 4 * vyi1 * vxi2 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi1 ^ 2 * vxi2 ^ 4 * x - 4 * vyi1 * vxi2 ^ 2 * vyi2 ^ 3 - 12 * vyi2 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 8 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi2 ^ 2 * vxi1 ^ 4 * x - 8 * vyi2 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi2 ^ 2 * vxi2 ^ 4 * x - 2 * vyi1 ^ 4 * vxi1 * vxi2 - vxi1 ^ 6 - vxi2 ^ 6 - 2 * vyi2 ^ 2 * vxi1 ^ 4 - vyi2 ^ 4 * vxi1 ^ 2 - 2 * vyi2 ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * vxi2 ^ 2 - 15 * vxi1 ^ 4 * vxi2 ^ 2 - 6 * vxi1 ^ 5 * vxi2 - 20 * vxi1 ^ 3 * vxi2 ^ 3 - 15 * vxi1 ^ 2 * vxi2 ^ 4 - 6 * vxi2 ^ 5 * vxi1 + 2 * vxi2 ^ 6 * x - vxi2 ^ 6 * x ^ 2 - 2 * vyi1 ^ 2 * vxi1 ^ 4 - vyi1 ^ 4 * vxi1 ^ 2 - 2 * vyi1 ^ 2 * vxi2 ^ 4 - vyi1 ^ 4 * vxi2 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi1 ^ 2 - 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f - 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f - 4 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 2 * vyi2 ^ 2 - 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vyi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * f ^ 2 * vyi2 ^ 2 + 8 * vyi1 * vxi1 ^ 3 * vxi2 * vyi2 * x + 8 * vyi1 * vxi1 * vxi2 ^ 3 * vyi2 * x + 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 * f + 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 * vxi1 ^ 2 * vyi2 * x * vxi2 ^ 2 + 8 * vyi1 ^ 2 * vxi1 * vxi2 * x * vyi2 ^ 2 + 8 * vyi1 ^ 3 * x * vyi2 * vxi1 * vxi2 + 8 * vyi2 ^ 3 * x * vyi1 * vxi1 * vxi2 - 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f * vyi1 ^ 2)) + (vyi1 ^ 3) + (vyi2 ^ 3) - (vyi1 ^ 3 * x) - (vyi2 ^ 3 * x) + (3 * vyi1 ^ 2 * vyi2) + (vyi1 * vxi1 ^ 2) + (vyi1 * vxi2 ^ 2) + (3 * vyi1 * vyi2 ^ 2) + (vyi2 * vxi1 ^ 2) + (vyi2 * vxi2 ^ 2)) / (vyi2 ^ 2 + vxi1 ^ 2 + 2 * vyi2 * vyi1 + vxi2 ^ 2 + vyi1 ^ 2 + 2 * vxi2 * vxi1) * vyi1 + ((2 * vyi1 * x * vxi1 ^ 2 * f) + (2 * vyi1 * x * vxi2 ^ 2 * f) + (2 * vyi1 * x * vyi2 ^ 2 * f) + (2 * vyi2 * x * vxi1 ^ 2 * f) + (2 * vyi2 * x * vxi2 ^ 2 * f) + (2 * vyi2 * x * vyi1 ^ 2 * f) + (2 * vyi2 ^ 3 * x * f) + (2 * vyi1 ^ 3 * x * f) + (2 * vyi1 * vxi1 * vxi2) - (vyi1 * x * vxi1 ^ 2) - (vyi1 * x * vxi2 ^ 2) - (vyi1 * x * vyi2 ^ 2) + (2 * vyi2 * vxi1 * vxi2) - (vyi2 * x * vxi1 ^ 2) - (vyi2 * x * vyi1 ^ 2) - (vyi2 * x * vxi2 ^ 2) + sqrt((-16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f ^ 2 + 16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f - 2 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 - 2 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 - 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f ^ 2 + 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f - 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f ^ 2 + 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi2 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi2 ^ 2 + 12 * vxi1 ^ 4 * x ^ 2 * f * vxi2 ^ 2 - 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 12 * vxi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi2 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 + 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f + 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 12 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi1 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 16 * vyi1 * vxi2 ^ 3 * vyi2 * vxi1 + 4 * vyi1 * vxi2 ^ 4 * vyi2 * x - 8 * vyi1 * vyi2 ^ 3 * vxi1 * vxi2 + 16 * vyi2 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi1 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi2 ^ 2 + 4 * vyi1 ^ 4 * x * vxi1 * vxi2 + 4 * vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 * f + 4 * vyi2 ^ 3 * x * vyi1 * vxi1 ^ 2 + 4 * vyi2 ^ 3 * x * vyi1 * vxi2 ^ 2 + 4 * vyi2 ^ 4 * x * vxi1 * vxi2 + 4 * vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 * f - 8 * vyi1 ^ 3 * vyi2 * vxi1 * vxi2 - 12 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 * vxi2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi1 ^ 2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi2 ^ 2 - 24 * vyi1 * vxi1 ^ 2 * vyi2 * vxi2 ^ 2 + 16 * vyi1 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 16 * vyi1 * vxi1 ^ 3 * vyi2 * vxi2 + 4 * vyi1 * vxi1 ^ 4 * vyi2 * x + 2 * vxi1 ^ 6 * x - vxi1 ^ 6 * x ^ 2 - 4 * vxi1 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi1 ^ 6 * x ^ 2 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 + 2 * vyi1 ^ 4 * x * vxi1 ^ 2 + 2 * vyi1 ^ 4 * x * vxi2 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 2 * vyi2 ^ 4 * x * vxi1 ^ 2 + 2 * vyi2 ^ 4 * x * vxi2 ^ 2 - vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 - 2 * vxi2 * vxi1 ^ 5 * x ^ 2 - 2 * vxi2 ^ 5 * vxi1 * x ^ 2 - 4 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 + 14 * vxi1 ^ 2 * vxi2 ^ 4 * x - 3 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 - 3 * vxi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 8 * vxi1 ^ 5 * vxi2 * x + 16 * vxi1 ^ 3 * vxi2 ^ 3 * x + 14 * vxi1 ^ 4 * x * vxi2 ^ 2 - 4 * vxi2 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi2 ^ 6 * x ^ 2 * f + 8 * vxi2 ^ 5 * vxi1 * x - 2 * vyi2 ^ 4 * vxi1 * vxi2 - 2 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi1 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi2 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi2 ^ 2 - 12 * vyi1 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 4 * vyi1 * vxi1 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi1 ^ 2 * vxi1 ^ 4 * x - 4 * vyi1 * vxi1 ^ 2 * vyi2 ^ 3 - 4 * vyi1 * vxi2 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi1 ^ 2 * vxi2 ^ 4 * x - 4 * vyi1 * vxi2 ^ 2 * vyi2 ^ 3 - 12 * vyi2 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 8 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi2 ^ 2 * vxi1 ^ 4 * x - 8 * vyi2 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi2 ^ 2 * vxi2 ^ 4 * x - 2 * vyi1 ^ 4 * vxi1 * vxi2 - vxi1 ^ 6 - vxi2 ^ 6 - 2 * vyi2 ^ 2 * vxi1 ^ 4 - vyi2 ^ 4 * vxi1 ^ 2 - 2 * vyi2 ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * vxi2 ^ 2 - 15 * vxi1 ^ 4 * vxi2 ^ 2 - 6 * vxi1 ^ 5 * vxi2 - 20 * vxi1 ^ 3 * vxi2 ^ 3 - 15 * vxi1 ^ 2 * vxi2 ^ 4 - 6 * vxi2 ^ 5 * vxi1 + 2 * vxi2 ^ 6 * x - vxi2 ^ 6 * x ^ 2 - 2 * vyi1 ^ 2 * vxi1 ^ 4 - vyi1 ^ 4 * vxi1 ^ 2 - 2 * vyi1 ^ 2 * vxi2 ^ 4 - vyi1 ^ 4 * vxi2 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi1 ^ 2 - 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f - 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f - 4 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 2 * vyi2 ^ 2 - 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vyi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * f ^ 2 * vyi2 ^ 2 + 8 * vyi1 * vxi1 ^ 3 * vxi2 * vyi2 * x + 8 * vyi1 * vxi1 * vxi2 ^ 3 * vyi2 * x + 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 * f + 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 * vxi1 ^ 2 * vyi2 * x * vxi2 ^ 2 + 8 * vyi1 ^ 2 * vxi1 * vxi2 * x * vyi2 ^ 2 + 8 * vyi1 ^ 3 * x * vyi2 * vxi1 * vxi2 + 8 * vyi2 ^ 3 * x * vyi1 * vxi1 * vxi2 - 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f * vyi1 ^ 2)) + (vyi1 ^ 3) + (vyi2 ^ 3) - (vyi1 ^ 3 * x) - (vyi2 ^ 3 * x) + (3 * vyi1 ^ 2 * vyi2) + (vyi1 * vxi1 ^ 2) + (vyi1 * vxi2 ^ 2) + (3 * vyi1 * vyi2 ^ 2) + (vyi2 * vxi1 ^ 2) + (vyi2 * vxi2 ^ 2)) / (vyi2 ^ 2 + vxi1 ^ 2 + 2 * vyi2 * vyi1 + vxi2 ^ 2 + vyi1 ^ 2 + 2 * vxi2 * vxi1) * vyi2 - (2 * vyi2 * vyi1) + (2 * vxi2 * vxi1) + (vxi1 ^ 2) - (vyi1 ^ 2) + (vxi2 ^ 2) - (vyi2 ^ 2) - (2 * x * vxi1 ^ 2 * f) - (2 * x * vxi2 ^ 2 * f) - (2 * x * vyi2 ^ 2 * f) - (2 * x * vyi1 ^ 2 * f) + (x * vxi1 ^ 2) + (vyi1 ^ 2 * x) + (vxi2 ^ 2 * x) + (vyi2 ^ 2 * x)) / (vxi1 + vxi2) / 0.2e1;
vyf2 = -((2 * vyi1 * x * vxi1 ^ 2 * f) + (2 * vyi1 * x * vxi2 ^ 2 * f) + (2 * vyi1 * x * vyi2 ^ 2 * f) + (2 * vyi2 * x * vxi1 ^ 2 * f) + (2 * vyi2 * x * vxi2 ^ 2 * f) + (2 * vyi2 * x * vyi1 ^ 2 * f) + (2 * vyi2 ^ 3 * x * f) + (2 * vyi1 ^ 3 * x * f) + (2 * vyi1 * vxi1 * vxi2) - (vyi1 * x * vxi1 ^ 2) - (vyi1 * x * vxi2 ^ 2) - (vyi1 * x * vyi2 ^ 2) + (2 * vyi2 * vxi1 * vxi2) - (vyi2 * x * vxi1 ^ 2) - (vyi2 * x * vyi1 ^ 2) - (vyi2 * x * vxi2 ^ 2) + sqrt((-16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f ^ 2 + 16 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 * f - 2 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 - 2 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 - 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f ^ 2 + 8 * vxi2 * vxi1 ^ 5 * x ^ 2 * f - 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f ^ 2 + 8 * vxi2 ^ 5 * vxi1 * x ^ 2 * f - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 * vxi1 ^ 3 * x ^ 2 * vyi2 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi1 ^ 2 - 4 * vxi2 ^ 3 * vxi1 * x ^ 2 * vyi2 ^ 2 + 12 * vxi1 ^ 4 * x ^ 2 * f * vxi2 ^ 2 - 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 12 * vxi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f ^ 2 + 12 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi2 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 + 8 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 * f + 8 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi1 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 12 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 * x + 12 * vyi1 ^ 2 * vxi1 * vxi2 ^ 3 * x - 4 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 * f - 16 * vyi1 * vxi2 ^ 3 * vyi2 * vxi1 + 4 * vyi1 * vxi2 ^ 4 * vyi2 * x - 8 * vyi1 * vyi2 ^ 3 * vxi1 * vxi2 + 16 * vyi2 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi1 ^ 2 - 4 * vyi2 ^ 4 * x ^ 2 * f ^ 2 * vxi2 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi1 ^ 2 + 4 * vyi1 ^ 3 * x * vyi2 * vxi2 ^ 2 + 4 * vyi1 ^ 4 * x * vxi1 * vxi2 + 4 * vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 * f + 4 * vyi2 ^ 3 * x * vyi1 * vxi1 ^ 2 + 4 * vyi2 ^ 3 * x * vyi1 * vxi2 ^ 2 + 4 * vyi2 ^ 4 * x * vxi1 * vxi2 + 4 * vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 * f + 4 * vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 * f - 8 * vyi1 ^ 3 * vyi2 * vxi1 * vxi2 - 12 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 * vxi2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi1 ^ 2 + 4 * vyi1 ^ 2 * vyi2 ^ 2 * x * vxi2 ^ 2 - 24 * vyi1 * vxi1 ^ 2 * vyi2 * vxi2 ^ 2 + 16 * vyi1 ^ 2 * vxi1 ^ 2 * x * vxi2 ^ 2 - 16 * vyi1 * vxi1 ^ 3 * vyi2 * vxi2 + 4 * vyi1 * vxi1 ^ 4 * vyi2 * x + 2 * vxi1 ^ 6 * x - vxi1 ^ 6 * x ^ 2 - 4 * vxi1 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi1 ^ 6 * x ^ 2 * f - 2 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 4 + 2 * vyi1 ^ 4 * x * vxi1 ^ 2 + 2 * vyi1 ^ 4 * x * vxi2 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi1 ^ 2 - vyi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 2 * vyi2 ^ 4 * x * vxi1 ^ 2 + 2 * vyi2 ^ 4 * x * vxi2 ^ 2 - vyi2 ^ 4 * x ^ 2 * vxi1 ^ 2 - 2 * vxi2 * vxi1 ^ 5 * x ^ 2 - 2 * vxi2 ^ 5 * vxi1 * x ^ 2 - 4 * vxi2 ^ 3 * vxi1 ^ 3 * x ^ 2 + 14 * vxi1 ^ 2 * vxi2 ^ 4 * x - 3 * vxi1 ^ 2 * x ^ 2 * vxi2 ^ 4 - 3 * vxi1 ^ 4 * x ^ 2 * vxi2 ^ 2 + 8 * vxi1 ^ 5 * vxi2 * x + 16 * vxi1 ^ 3 * vxi2 ^ 3 * x + 14 * vxi1 ^ 4 * x * vxi2 ^ 2 - 4 * vxi2 ^ 6 * x ^ 2 * f ^ 2 + 4 * vxi2 ^ 6 * x ^ 2 * f + 8 * vxi2 ^ 5 * vxi1 * x - 2 * vyi2 ^ 4 * vxi1 * vxi2 - 2 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 4 - 2 * vyi2 ^ 2 * x ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * x ^ 2 * vxi2 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi1 ^ 2 - 4 * vyi1 ^ 3 * vyi2 * vxi2 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi1 ^ 2 - 6 * vyi1 ^ 2 * vyi2 ^ 2 * vxi2 ^ 2 - 12 * vyi1 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 4 * vyi1 * vxi1 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi1 ^ 2 * vxi1 ^ 4 * x - 4 * vyi1 * vxi1 ^ 2 * vyi2 ^ 3 - 4 * vyi1 * vxi2 ^ 4 * vyi2 - 8 * vyi1 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi1 ^ 2 * vxi2 ^ 4 * x - 4 * vyi1 * vxi2 ^ 2 * vyi2 ^ 3 - 12 * vyi2 ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 - 8 * vyi2 ^ 2 * vxi1 ^ 3 * vxi2 + 4 * vyi2 ^ 2 * vxi1 ^ 4 * x - 8 * vyi2 ^ 2 * vxi2 ^ 3 * vxi1 + 4 * vyi2 ^ 2 * vxi2 ^ 4 * x - 2 * vyi1 ^ 4 * vxi1 * vxi2 - vxi1 ^ 6 - vxi2 ^ 6 - 2 * vyi2 ^ 2 * vxi1 ^ 4 - vyi2 ^ 4 * vxi1 ^ 2 - 2 * vyi2 ^ 2 * vxi2 ^ 4 - vyi2 ^ 4 * vxi2 ^ 2 - 15 * vxi1 ^ 4 * vxi2 ^ 2 - 6 * vxi1 ^ 5 * vxi2 - 20 * vxi1 ^ 3 * vxi2 ^ 3 - 15 * vxi1 ^ 2 * vxi2 ^ 4 - 6 * vxi2 ^ 5 * vxi1 + 2 * vxi2 ^ 6 * x - vxi2 ^ 6 * x ^ 2 - 2 * vyi1 ^ 2 * vxi1 ^ 4 - vyi1 ^ 4 * vxi1 ^ 2 - 2 * vyi1 ^ 2 * vxi2 ^ 4 - vyi1 ^ 4 * vxi2 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi1 ^ 2 + 16 * vxi2 * vxi1 ^ 3 * x ^ 2 * f * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi2 ^ 2 - 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 ^ 3 * vxi1 * x ^ 2 * f * vyi1 ^ 2 - 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 4 * f - 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f ^ 2 + 8 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 4 * f - 4 * vxi2 * vxi1 * x ^ 2 * vyi1 ^ 2 * vyi2 ^ 2 - 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vxi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * f ^ 2 * vyi2 ^ 2 - 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * f ^ 2 * vyi2 ^ 2 + 8 * vyi1 * vxi1 ^ 3 * vxi2 * vyi2 * x + 8 * vyi1 * vxi1 * vxi2 ^ 3 * vyi2 * x + 16 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi1 ^ 2 * vyi2 ^ 2 * f + 8 * vyi1 ^ 2 * x ^ 2 * vxi2 ^ 2 * vyi2 ^ 2 * f + 16 * vyi2 ^ 2 * x ^ 2 * vxi1 ^ 2 * vxi2 ^ 2 * f + 8 * vyi1 * vxi1 ^ 2 * vyi2 * x * vxi2 ^ 2 + 8 * vyi1 ^ 2 * vxi1 * vxi2 * x * vyi2 ^ 2 + 8 * vyi1 ^ 3 * x * vyi2 * vxi1 * vxi2 + 8 * vyi2 ^ 3 * x * vyi1 * vxi1 * vxi2 - 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f ^ 2 * vyi1 ^ 2 + 16 * vxi2 * vxi1 * x ^ 2 * vyi2 ^ 2 * f * vyi1 ^ 2)) + (vyi1 ^ 3) + (vyi2 ^ 3) - (vyi1 ^ 3 * x) - (vyi2 ^ 3 * x) + (3 * vyi1 ^ 2 * vyi2) + (vyi1 * vxi1 ^ 2) + (vyi1 * vxi2 ^ 2) + (3 * vyi1 * vyi2 ^ 2) + (vyi2 * vxi1 ^ 2) + (vyi2 * vxi2 ^ 2)) / (vyi2 ^ 2 + vxi1 ^ 2 + 2 * vyi2 * vyi1 + vxi2 ^ 2 + vyi1 ^ 2 + 2 * vxi2 * vxi1) / 0.2e1 + vyi1 + vyi2;

if ~(isreal(vyf1) && isreal(vyf2) && isreal(vxf1) && isreal(vxf2))
    fprintf('Problem!  Fractionation not achievable.\n');
    fail = 1;
end

%checkf = (vxf1 ^ 2 + vyf1 ^ 2) / (vxf1 ^ 2 + vyf1 ^ 2 + vxf2 ^ 2 + vyf2 ^ 2);
%fprintf(1, 'Final energy fraction in body 1 is %f.\n', checkf);


