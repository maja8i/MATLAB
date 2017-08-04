%Generic code to clear all variables and stop writing to diary file. Use in
%between different runs of Bezier_volumes_main.m, when the previous code
%did not run to the end (e.g., was stopped during debugging).

diary off
fclose all;
close all
clear
clc