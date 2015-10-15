% Title : Computes the PSD matrix
% File  : compute_psd_matrix.m
% -------------------------------------------------------------------------
% Description :
% Computes the oversampled power spectral density
% of the amplified transmit signal of a MIMO system that uses maximum-ratio
% precoding over a frequency-selective channel specified by CHANNEL by pulse
% shape filtering with PULSE.  OVERSAMPLING_FACTOR specifies the time-domain
% sampling rate i units of the baud-rate.
% ------------------------------------------------------------------------- 
% Revisions   :
%   Date       Version  Author  Description
%   15-Oct-15  1.0      chrism  created the file
% -------------------------------------------------------------------------                            
%   Author: Christopher Mollen   
% ------------------------------------------------------------------------- 

function [ psd ] = compute_psd_matrix( channel, pulse, oversamp_factor )

[nr_users, nr_antennae, nr_taps] = size(channel);

b_1 = 1;
b_2 = (-0.0373 + 1i*0.0065) / (1.0707 - 1i*0.0129);
b_1s = b_1 * ones(nr_antennae, 1);
b_2s = b_2 * ones(nr_antennae, 1);

backoff_dB = 0;
backoff = 10^(backoff_dB/10);

baudrate = 1;

agg_pulse = falta(pulse, pulse);

[ds_channel, nr_noncausal_samps] = downsample_channel(channel, pulse, oversamp_factor);
ds_channel = ds_channel * sqrt(oversamp_factor / sum(abs(agg_pulse).^2));
nr_ds_taps = length(ds_channel(1,1,:));

precoding_mat = zeros(nr_antennae, nr_users, nr_ds_taps);
for tap_id = 1:nr_ds_taps
    % maximum-ratio precoding
    precoding_mat(:, :, tap_id) = ds_channel(:,:, nr_ds_taps - tap_id + 1)' / sqrt(nr_antennae);
end

disc_time_preamp_corr = zeros(nr_antennae, nr_antennae, nr_ds_taps);
for shift = 0:(nr_ds_taps - 1)
    for tap_id = 1:(nr_ds_taps - shift)
        disc_time_preamp_corr(:, :, shift + 1) = disc_time_preamp_corr(:, :, shift + 1) ...
            + conj(precoding_mat(:, :, tap_id)) * precoding_mat(:, :, shift + tap_id).' / nr_users;
    end
end

agg_pulse_zero_id = ceil(length(agg_pulse) / 2);

nr_nonneg_samples = 4 * nr_taps; %nr of samples computed for preamp_corr (only positive samples computed, because of symmetry)
preamp_corr = zeros(nr_antennae, nr_antennae, nr_nonneg_samples);

for time_id = 1:nr_nonneg_samples
    for time = (-nr_ds_taps+1):(nr_ds_taps-1)
        if agg_pulse_zero_id + time_id - 1 - time * oversamp_factor <= length(agg_pulse)
            if time < 0
                preamp_corr(:, :, time_id) = preamp_corr(:, :, time_id) ...
                    + baudrate * disc_time_preamp_corr(:, :, -time + 1)' * agg_pulse(agg_pulse_zero_id + time_id - 1 - time * oversamp_factor);
            elseif time == 0
                preamp_corr(:, :, time_id) = preamp_corr(:, :, time_id) ...
                    + baudrate * disc_time_preamp_corr(:, :, 1) * agg_pulse(agg_pulse_zero_id + time_id - 1);
            elseif time > 0
                preamp_corr(:, :, time_id) = preamp_corr(:, :, time_id) ...
                    + baudrate * disc_time_preamp_corr(:, :, time + 1) * agg_pulse(agg_pulse_zero_id + time_id - 1 - time * oversamp_factor);
            end
        end
    end
end

av_power = mean(diag(preamp_corr(:, :, 1)));
preamp_corr = preamp_corr / av_power; % unity per-antenna power before amplifier.

% find the largest positive one-dB compression point
onedB_comp_1 = (-real(b_2 / b_1) + sqrt((real(b_2 / b_1))^2 - (1 - 10^(-.1)) * abs(b_2 / b_1)^2)) / abs(b_2 / b_1)^2;
onedB_comp_2 = (-real(b_2 / b_1) - sqrt((real(b_2 / b_1))^2 - (1 - 10^(-.1)) * abs(b_2 / b_1)^2)) / abs(b_2 / b_1)^2;
if onedB_comp_2 <= 0
    onedB_comp = onedB_comp_1;
else
    onedB_comp = onedB_comp_2;
end

op_point = onedB_comp * backoff;
preamp_corr = preamp_corr * op_point; % set the operating point specified by the backoff.

amp_corr = zeros(nr_antennae, nr_antennae, nr_nonneg_samples);
for antenna1_id = 1:nr_antennae
    for antenna2_id = 1:nr_antennae
        for sample_id = 1:nr_nonneg_samples
            amp_corr(antenna1_id, antenna2_id, sample_id) = ...
                conj(b_1s(antenna1_id)) * b_1s(antenna2_id) * preamp_corr(antenna1_id, antenna2_id, sample_id)...
                + 2 * preamp_corr(antenna1_id, antenna2_id, sample_id) ...
                * ( conj(b_1s(antenna1_id)) * b_2s(antenna2_id) * preamp_corr(antenna1_id, antenna1_id, 1) ...
                + conj(b_2s(antenna1_id)) * b_1s(antenna2_id) * preamp_corr(antenna2_id, antenna2_id, 1) ...
                + conj(b_2s(antenna1_id)) * b_2s(antenna2_id) ...
                * ( 2 * preamp_corr(antenna1_id, antenna1_id, 1) * preamp_corr(antenna2_id, antenna2_id, 1) ...
                + abs(preamp_corr(antenna1_id, antenna2_id, sample_id))^2));
        end
    end
end

nr_freqs = 2 * nr_nonneg_samples - 1; % the chosen freq resolution
psd = zeros(nr_antennae, nr_antennae, nr_freqs);

for freq_id = 1:nr_freqs
    freq = freq_id - 1;
    psd(:, :, freq_id) = amp_corr(:, :, 1);
    
    for sample_id = 2:nr_nonneg_samples
        time = sample_id - 1;
        psd(:, :, freq_id) = psd(:, :, freq_id) ...
            + amp_corr(:, :, sample_id) * exp(-2j * pi * time * freq / nr_freqs) + amp_corr(:, :, sample_id)' * exp(2j * pi * time * freq / nr_freqs);
    end
end
end

