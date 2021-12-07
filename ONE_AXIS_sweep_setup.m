%% Run Sweeps / Acquire data from lock-in amplifier 
%amplitude1 = 500e-3;  %drive amplitude (V)
%amplitude2 = 500e-3;
close all
clear
clc
%{
sample_rate = 225;%Hz
demod_filt_order = 4;
desired_demod_BW = 10.2 Hz; 
demod_TC = compute_demod_TC(demod_filt_order,desired_demod_BW); %refer to
page 776 of old HF2LI manual
%write function that converts order and desired 3db BW to Time constant
% 
%}
%axis2 scale up = 2.05276(off-resonance)or 4.667 (on resonance time domain)
amplitude = [0.1 (2.05*0.1)];%500e-3;  %drive amplitude (V) (MAXIMUM 1 V amplitude!)  
num_points = 1024;%1024; %1024
TC_settle = 25;
settle_time = 0.75;%0.75 :0.1(2.5min,512pts);0.5(9.5 min,512pts);1(16min,512pts)
midpoint_freq = 10538.980450;% (10541.1) 10553;
freq_span = 10; %fine sweep = 2,coarse sweep = 20  %freq span will be +/- this value 
start_f = midpoint_freq - freq_span;%10530;  %(wide sweep 10490 - 10540) (narrow sweep: 10520 - 10560)
stop_f =  midpoint_freq + freq_span;%10570;
samp_avg = 8;

filter_slope = 24; %dB/Oct  6(1st),12(2nd),18(3rd),24(4th),30(5th),36(6th),42(7th),48(8th)
filter_BW = 10.2;%1.95;%1.95; %Hz (1.95,4.86,10.2, 
demod_samp_rate = 225; %Hz (225,450, ...)

scan_order = 0; %sequential = 0, binary = 1, bidirectional (forward&reverse) = 2, reverse = 3
forward_and_reverse = 0;
sweep_which_axis = 2; % drive 180 (sense0/45) = 0 , drive 225 (sense0/45) = 1, drive 180 and drive 225 = 2
diff_sense = 0; % 0 = single ended sense, 1 = differential sensing on
sense_range = 0.5; %sensing input range (vpk) 
aux1_out = 0; %G3.1 (12) 22.5 deg axis
aux2_out = 0; %G2 (4) 45 deg axis
aux3_out = 0; %G4.2 (8) -22.5 deg axis
%====================================================
base_file_loc = 'F:\Users\cwboyd\Desktop\DATA_2021\May_Sept_21\gain_phase_logs\';
logfile_name = '7_13_21.txt';
file_name = strcat(base_file_loc,logfile_name);
settings_filename = 'base_sweep_settings_7_22_19.xml';
step_num = 2;
pause_time = 20; %seconds
%DRIVE: 202.5/157.5 (180/225 = 0V)
%DRIVE: 180/225 (202.5(G3.3)/157.5(G4.4) = 0V)
comment1 = 'DRIVE: 180/225 (202.5(G3.3(9))/157.5(G4.4(6)) = 0V) Device# BRG_12: BIAS = +150V, G1 = -155.25, G3.4 = 0, G2 = 0V';
comment2 =  'G3.1 = G4.2 = 0V, G4.4 = G3.3 = 0V, rest = +150 V, metal down - SOUTH, ASIC TIA 0x80808080';%'Discrete TIA gain = 0.5 MOhm?';
%G1 = -134.865V, G3.4 = +45.190V, G2 = 0V, G3.1 = 0V, G4.2 = 0V, Bias = +200V, rest = +200V; 180 drive2P(20) 225 drive1P(23)
%comment1 = 'bias=+200V,Tchuck=233.15,Tcontrol=ON, official measurement';
%comment2 = ' 180 drive2P(20) and drive1P(23) OR 202.5 drive drive2P(9) and drive1P(6) ';