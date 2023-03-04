clear; close all; clc;
%File for loading shot data from hdf5 files
time = h5read("220926010.hdf5","/time");
%Flux coil currents for a plasma
i_fcoil_1 = h5read("220926010.hdf5","/inj/i1inj");
i_fcoil_2 = h5read("220926010.hdf5","/inj/i2inj");
i_fcoil_3 = h5read("220926010.hdf5","/inj/i3inj");
i_fcoil_4 = h5read("220926010.hdf5","/inj/i4inj");

%i_tor

i_tor = h5read("220926010.hdf5","/itor");


