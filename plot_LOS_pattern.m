clear;

nr_users = 10;
nr_antennae = 100;
oversamp_factor = 5;
nr_taps = oversamp_factor * 15;
pulse_type = 'RRC';
rolloff = 0.22;
pulse = rcosine(1, oversamp_factor, 'sqrt', rolloff, 300).'/ sqrt(oversamp_factor);
fading_model = 'LOS';
antenna_spacing_wl = 1/2;

% Generate channel realisation
incidence_angles = rand(nr_users,1) * pi - pi / 2;
if strcmp(fading_model, 'Rayleigh')
    channel = sqrt(0.5 / nr_taps) * randn(nr_users, nr_antennae, nr_taps) + 1j * sqrt(0.5 / nr_taps) * randn(nr_users, nr_antennae, nr_taps);
elseif strcmp(fading_model, 'LOS')
    channel = zeros(nr_users, nr_antennae, nr_taps);
    for user_id = 1:nr_users
        incidence_angle = incidence_angles(user_id);
        phase_rot = rand(1) * 2 * pi;
        stering_vec = exp(2j * pi * ((1:nr_antennae) - 1)' * antenna_spacing_wl * sin(incidence_angle));
        channel(user_id, :, 1) = exp(1j * phase_rot) * stering_vec;
    end
end

psd = compute_psd_matrix(channel, pulse, oversamp_factor);

% Compute dominant eigenvalues at different frequencies
nr_freqs = length(psd(1, 1, :));
freqs = linspace(-oversamp_factor / 2, oversamp_factor / 2, nr_freqs);
dom_eigs = zeros(nr_freqs, 1);

for freq_id = 1:nr_freqs
    dom_eigs(freq_id) = nr_antennae * eigs(psd(:, :, freq_id), 1);
end

% Compute Pobmax
bandwidth_samp = (1 + rolloff) * nr_freqs / oversamp_factor;
rightband = [(floor(bandwidth_samp / 2) + 1):(floor(3 * bandwidth_samp / 2))];
leftband = [(nr_freqs - 3 * floor(bandwidth_samp / 2) + 1):(nr_freqs - floor(bandwidth_samp / 2))];
Pobmax = max(sum(real(dom_eigs(leftband))), sum(real(dom_eigs(rightband))));

% Compute LOS pattern
angles = linspace(-pi/2, pi/2, 500);
nr_angles = length(angles);
oob_powers = zeros(nr_angles, nr_freqs);
for angle_id = 1:nr_angles
    for freq_id = 1:nr_freqs        
        ang = angles(angle_id);
        channel_vec = exp(2j * pi * ((1:nr_antennae) - 1)' * sin(ang) * antenna_spacing_wl);
        oob_powers(angle_id, freq_id) = channel_vec' * psd(:, :, freq_id) * channel_vec;
    end
end
Pob = max(sum(real(oob_powers(:, leftband)), 2), sum(real(oob_powers(:, rightband)), 2));

% Plot
figure(4)

antenna_id = 1;
tx_psd = zeros(nr_freqs, 1);
for freq_id = 1:nr_freqs
    tx_psd(freq_id) = trace(psd(:,:,freq_id));
end
tx_Pob = max(sum(real(tx_psd(leftband))), sum(real(tx_psd(rightband))));

ymax = 10 * log10(Pobmax) + 1.5;
ymin = min(10 * log10(Pob)) - 1;

plot(angles, 10 * log10(Pob), 'LineWidth', 1)
hold on
lines_x = [incidence_angles, incidence_angles]';
lines_y = [ones(1, nr_users) * ymin; ones(1, nr_users) * ymax];
plot(lines_x, lines_y, 'r')
plot([-pi/2, pi/2], 10 * log10([Pobmax Pobmax]));
text(-pi/2 + .1, 10 * log10(Pobmax) - 0.3, 'P_{ob,max}', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
h = plot([-pi/2, pi/2], 10 * log10([tx_Pob tx_Pob]), '--');
hold off

xlim([-pi/2, pi/2])
ylim([ymin, ymax])
xlabel('Incidence Angle \theta [rad]', 'FontSize', 15)
ylabel('P_{ob}(\theta) [dB]', 'FontSize', 15)
legend(h,'radiated out-of-band power', 'Location', 'SE')
