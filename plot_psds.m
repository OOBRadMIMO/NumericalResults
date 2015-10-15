clear

nr_users = 10;
nr_antennae = 100;
oversamp_factor = 5;
nr_taps = oversamp_factor * 15;
pulse_type = 'RRC';
rolloff = 0.22;
pulse = rcosine(1, oversamp_factor, 'sqrt', rolloff, 300).' /  sqrt(oversamp_factor);
fading_model = 'Rayleigh';
antenna_spacing_wl = 1/2;

% Generate a channel realisation
incidence_angles = rand(nr_users,1) * pi - pi / 2;
if strcmp(fading_model, 'Rayleigh')
    channel = sqrt(0.5 / nr_taps) * randn(nr_users, nr_antennae, nr_taps) + 1j * sqrt(0.5 / nr_taps) * randn(nr_users, nr_antennae, nr_taps);
elseif strcmp(fading_model, 'LOS')
    channel = zeros(nr_users, nr_antennae, nr_taps);
    for user_id = 1:nr_users
        incidence_angle = incidence_angles(user_id);
        phase_rot = rand(1) * 2 * pi;
        stering_vec = exp(2j * pi * ((1:nr_antennae)-1)' * antenna_spacing_wl * sin(incidence_angle));
        channel(user_id, :, 1) = exp(1j * phase_rot) * stering_vec;
    end
end

psd = compute_psd_matrix(channel, pulse, oversamp_factor);

% Compute the dominant eigenvalue at all frequencies
nr_freqs = length(psd(1, 1, :));
freqs = linspace(-oversamp_factor / 2, oversamp_factor / 2, nr_freqs);
dom_eigs = zeros(nr_freqs, 1);
for freq_id = 1:nr_freqs
    dom_eigs(freq_id) = nr_antennae * eigs(psd(:,:,freq_id),1);
end

% Generate random channel for a victim, independent of user channels in channel
if strcmp(fading_model, 'Rayleigh')
    rand_user_channel = sqrt(0.5 / nr_taps) * randn(nr_antennae, nr_taps) + 1j * sqrt(0.5 / nr_taps) * randn(nr_antennae, nr_taps);
elseif strcmp(fading_model, 'LOS')
     phase_rot = rand(1) * 2 * pi;
     incidence_angle = rand(1, 1) * pi - pi/2;
     stering_vec = exp(2j * pi * ((1:nr_antennae) - 1)' * antenna_spacing_wl * sin(incidence_angle));
     rand_user_channel = exp(1j*phase_rot) * stering_vec;
end

% Compute channels in frequency domain
rand_user_channel_freq = fft(rand_user_channel, nr_freqs, 2);
freq_channel = fft(channel, nr_freqs, 3);

% Compute the received PSDs of the users and the victim
bandwidth_samp = (1 + rolloff) * nr_freqs / oversamp_factor;
inband = [1:floor(bandwidth_samp / 2), (nr_freqs - floor(bandwidth_samp / 2) + 1):nr_freqs];
user_PSDs = zeros(nr_freqs, nr_users);
rand_user_PSD = zeros(nr_freqs, 1);
Pib = zeros(nr_users, 1);
for user_id = 1:nr_users
    for freq_id = 1:nr_freqs
        h = freq_channel(user_id, :, freq_id);
        h = h(:);
        user_PSDs(freq_id, user_id) = h' * psd(:, :, freq_id) * h;
        h = rand_user_channel_freq(:, freq_id);
        rand_user_PSD(freq_id) = h' * psd(:,:,freq_id) * h;
    end
    Pib(user_id) = sum(user_PSDs(inband, user_id));
end
[Pibmin, min_user_id]  = min(real(Pib));

% Plot
figure(4)

y_vals = fftshift(10 * log10(max(0, real(dom_eigs))));
sa_psd = zeros(nr_freqs, 1);
for antenna_id = 1:nr_antennae
    psd_one_antenna = max(0,real(psd(antenna_id,antenna_id,:)));
    sa_psd = sa_psd + psd_one_antenna(:);
end

plot(freqs, y_vals)
hold on
user_id = min_user_id;
plot(freqs, fftshift(10 * log10(user_PSDs(:, user_id))))
plot(freqs, fftshift(10 * log10(rand_user_PSD)))
plot(freqs, fftshift(10 * log10(sa_psd(:))))

xlabel('Normalized Frequency fT')
ylabel('Power [dB]')
lower_lim = -10;
upper_lim = 45;
ylim([lower_lim, upper_lim])
plot([(-0.5 - rolloff/2) (0.5 + rolloff/2)], [lower_lim lower_lim], 'LineWidth', 3)
legend({'maximum PSD', 'user with min $P_\text{ib}(\mathbf{\theta}_k)$', 'random spatial point', 'sum per-antenna PSD'}, 'Position', [0.01, 0.08, 0.25, 0.2])
text(0, lower_lim, 'in-band', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

hold off

