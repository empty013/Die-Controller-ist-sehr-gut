clear all; close all
addpath("controller\")
global controller_func
controller_func = @u_2inp_p;
simulate_controlled_singletrack(10)