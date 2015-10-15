clear

nr_users = 1;
nr_antennae = 100;
oversamp_factor = 5;
nr_taps = oversamp_factor * 15;
rolloff = 0.22;
pulse = rcosine(1, oversamp_factor, 'sqrt', rolloff, 300).'/ sqrt(oversamp_factor);

fading_model = 'Rayleigh';

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
        stering_vec = exp(2j * pi * ((1:nr_antennae)-1)' * antenna_spacing_wl * sin(incidence_angle));
        channel(user_id, :, 1) = exp(1j * phase_rot) * stering_vec;
    end
end

psd = compute_psd_matrix(channel, pulse, oversamp_factor);
nr_freqs = length(psd(1, 1, :));

% Plot
figure(34)

bandwidth = 1 + rolloff;
studied_freqs = [bandwidth 3*bandwidth/4 bandwidth/2 bandwidth/4 0];
studied_freq_ids = round(studied_freqs * nr_freqs / oversamp_factor) + 1;
nr_studied_freqs = length(studied_freqs);

eigenvalues = zeros(nr_antennae, nr_studied_freqs);
colors = {[.54, .5, .82], [0, .81, .71], [.41, .49, .57], [1, .39, .26], [0, .72, .90]};

leg= {};
for studied_freq_id = 1:nr_studied_freqs
    freq_id = studied_freq_ids(studied_freq_id);
    eigenvalues(:, studied_freq_id) = real(eigs(psd(:, :, freq_id), nr_antennae));
    
    [nr_occurences, eigenvalue_bins] = hist(eigenvalues(:, studied_freq_id), 100 * nr_antennae);
    probabilities = 1 - cumsum(nr_occurences) / nr_antennae;
    
    x_vals = 10 * log10(max(0, eigenvalue_bins));
    y_vals = probabilities;
    
    % Remove NaN and Inf
    x_vals_clean = x_vals(y_vals ~= NaN);
    y_vals_clean = y_vals(y_vals ~= NaN);
    y_vals_clean = y_vals_clean(x_vals_clean ~= -Inf);
    x_vals_clean = x_vals_clean(x_vals_clean ~= -Inf);
        
    hold on
    
    h(studied_freq_id) = plot(x_vals_clean, y_vals_clean, 'Color', colors{studied_freq_id});
    leg{studied_freq_id} = ['fT = ' num2str(round((freq_id - 1) * oversamp_factor / nr_freqs * 100) / 100)];
    
    mean_power_dB = 10 * log10(real(trace(psd(:, :, freq_id))) / nr_antennae);
    mean_prob = interp1(x_vals_clean, y_vals_clean, mean_power_dB);
    
    plot(mean_power_dB, mean_prob, '.', 'Color', colors{studied_freq_id}, 'MarkerSize', 25)
    
end

legend(h, leg)
xlabel('Eigenvalues of S_{yy}(f) [dB]')
ylabel('Fraction of Eigenvalues')

