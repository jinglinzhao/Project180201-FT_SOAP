% Modify from Project180131-FT_SOAP in order to pick up radnom phases
% instead of consecutive observations

% Use simulated spectra with planets AND stellar jitter: 
% /Volumes/DataSSD/SOAP_2/outputs/02.01/

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
SN              = 10000;
N_FILE          = 200;                               % number of CCF files
grid_size       = 0.1;
v0              = (-20 : grid_size : 20)';          % km/s
dir1            = '/Volumes/DataSSD/SOAP_2/outputs/02.01/';
dir2            = '/Volumes/DataSSD/SOAP_2/outputs/02.01/CCF_dat/';
RV              = importdata([dir1, 'RV.dat']) / 1000;      % activity induced RV [km/s]
idx             = (v0 > -10) & (v0 < 10);
v1              = v0(idx);
period_min      = 2;
period_max      = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%

% estimate the size of array FFT_power
filename    = [dir2, 'CCF', num2str(1), '.dat'];
A           = 1 - importdata(filename);
A           = A(idx);
[aa, bb]    = FUNCTION_FFT_noise(A, 0.1);
size1       = length(bb);
FFT_power   = zeros(size1, N_FILE);

% randomly select 200 different integers from 0 to 999
% 100 corresponds to one solar roation period ~ 25 days 
MJD = sort(randperm(1000,200))-1;

for n = 1:200
        
    i           = mod(MJD(n), 100);
    v_planet    = 10 * sin(i/100*0.618*2*pi + 1) * 0.001;         % km/s
    filename    = [dir2, 'CCF', num2str(i), '.dat'];
    A           = 1 - importdata(filename);
    A_spline    = spline(v0, A, v1+v_planet);
    
    [FFT_frequency, FFT_power(:, n)] = FUNCTION_FFT_noise(A_spline, 0.1);

    if 0            % Plot line profile in frequency space
        h_fft = figure;
        hold on 
        plot(FFT_frequency(1:int16(size1/3)), FFT_power(1:int16(size1/3), n))
        plot(FFT_frequency(1:int16(size1/3)), FFT_power(1:int16(size1/3), n), '.', 'markers',12)
        hold off
        xlabel('FT frequency domain')
        ylabel('Normalized Power')        
        title_name = ['FT - file', num2str(n)];
        out_eps = [title_name, '.eps'];
        print(out_eps, '-depsc')
        close(h_fft); 
    end
    
end     

% estimate the size of array Z_fft
[pxx_fft,f_fft] = plomb(FFT_power(1, :)', MJD/100*25, 0.5,  'normalized');
Z_fft = zeros(length(FFT_frequency), length(f_fft));

idx = FFT_power(:,1) > 0.01;
idx_array = 1:size1;
idx_max = max(idx_array(idx));

for j = 1:idx_max
    
    [pxx_fft,f_fft] = plomb(FFT_power(j, :)', MJD/100*25, 0.5,  'normalized');
    Z_fft(j,:) = pxx_fft;
    [pks_fft,locs_fft] = findpeaks(pxx_fft, f_fft);                         % find all the peaks in (pxx, f)
    [pks_maxs_fft, idx_maxs_fft] = sort(pks_fft, 'descend');                % sort "pks" in descending order; mark the indexies     
        
    h = figure;
        semilogx(1./f_fft, pxx_fft)
        xlim([period_min period_max])
        
        hold on
            % label the first 10 strongest signals
            for i = 1:10
                % Mark coefficient peaks (red)
                x = locs_fft(idx_maxs_fft(i));                                          % locations of the largest peaks -> harmonics
                y = pks_maxs_fft(i);
                if (1/x < period_max) && (1/x > period_min)
                    text(1/x, y, ['\leftarrow', num2str(1/x, '%3.2f')], 'fontsize', 8);
                end
            end

            xlabel('Period')
            ylabel('Normalized Power')    
            title_name = ['Frequency number', num2str(FFT_frequency(j))];
            title(title_name);
        hold off    
        
        out_eps = [title_name, '.eps'];
        print(out_eps, '-depsc')
    close(h);    
    
end        

x_idx =  (f_fft < 1/period_min);
[X_fft,Y_fft] = meshgrid(f_fft(x_idx), FFT_frequency(idx));
h = pcolor(X_fft,Y_fft, (Z_fft(idx,x_idx)));
set(h, 'edgecolor', 'none');
colorbar;
xlabel('Activity Frequency [1/days]');
ylabel('FT frequency domain');
title('Stacked periodogram');
print('Stacked periodogram', '-depsc')