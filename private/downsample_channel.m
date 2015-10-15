% Title : Downsample of the channel
% File  : cownsample_channel.m
% -------------------------------------------------------------------------
% Description :
% This script performs down-sampling of the channel.
% ------------------------------------------------------------------------- 
% Revisions   :
%   Date       Version  Author  Description
%   15-Oct-15  1.0      chrism  created the file
% -------------------------------------------------------------------------                            
%   Author: Christopher Mollen   
% ------------------------------------------------------------------------- 

function [ds_channel, delay] = downsample_channel(channel, pulse, oversamp_factor)
[nr_users, nr_antennae, nr_taps] = size(channel);
nr_samples = 2 * length(pulse) - 1;
delay = 2;
start_tap = -delay;
stop_tap = ceil(nr_taps / oversamp_factor) + delay - 1;

ds_channel = zeros(nr_users, nr_antennae, stop_tap - start_tap + 1);
for user_id = 1:nr_users
    for antenna_id = 1:nr_antennae
        sig = channel(user_id, antenna_id, :);
        sig = falta(pulse, falta(pulse, sig(:)));
        start_id = floor(nr_samples/2) + 1 + oversamp_factor * start_tap;
        stop_id = start_id + oversamp_factor * (stop_tap - start_tap);
        ds_channel(user_id, antenna_id, :) = sig(start_id:oversamp_factor:stop_id);
    end
end


end