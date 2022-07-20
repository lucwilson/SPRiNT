%% Code to simulate n (10) time series with the same spectrum as simulation challenge 1
len = 60;                               % len: simulation length (s)
fs = 200;                               % fs: sampling rate
ap_init = 1.5;                          % ap_init: initial aperiodic exponent (Hz^-1)
pk_cf = [8 18 8 8];                     % pk_cf: initial peak center frequencies (Hz)    
pk_amp = [1.2 0.9 1.2 1.2];             % pk_amp: initial peak amplitudes (a.u.)  
pk_sd = [1.2 1.4 1.2 1.2];              % pk_sd: initial peak standard deviations (Hz)  
pk_times = [8 32;15 10;41 5;47 5];      % pk_times: peak [start_time duration] (s)
pk_cf_chng = [0 0.3 0.8;-3 0.3 0.7;...  % pk_cf_chng: peak [change_in_cf change_start change_end] 
    0 0.3 0.8;0 0.3 0.8];               % in [(Hz) (prop_of_peak_time prop_of_peak_time)] 
pk_sd_chng = [0 0.3 0.8;0 0.3 0.8;...   % pk_sd_chng: peak [change_in_sd change_start change_end] 
    0 0.3 0.8;0 0.3 0.8];               % in [(Hz) (prop_of_peak_time prop_of_peak_time)] (unused)
ap_chng = [0.5 0.40 0.60];              % ap_chng: aperiodic [change_in_exponent change_start change_end] 
                                        % in [(Hz^-1) (prop_of_simulation prop_of_simulation)]
rot_f = 200;                            % rot_f: rotation frequency (Hz)
n = 10;                                 % n: number of simulations

xs1 = sim_combined_mat5(len, fs, ap_init, pk_cf, pk_amp, pk_sd, pk_times,...
    pk_cf_chng, pk_sd_chng, ap_chng, rot_f, n);


%% Code to simulate n (10) time series using randomly generated spectral parameters (new challenge 2 data)

len = 60;                               % len: simulation length (s)
ap_init_ranges = [0.8 2.2; 0.1 100];    % ap_init_ranges: initial aperiodic range
                                        % First row: aperiodic exponent range (Hz^-1)
                                        % Second row: rotation frequency range (Hz)
ap_chng_range = [-0.5 0.5];             % ap_chng_range: change range of aperiodic exponent (Hz^-1)
ap_chng_start_range = [0.2 0.6];        % ap_chng_start_range: start point of aperiodic shift, 
                                        % in proportion of recording (e.g., 60s recording: 0.2 = 12 s)
ap_chng_dur_range = [0.1 0.4];          % ap_chng_dur_range: duration range of aperiodic shift, 
                                        % in proportion of recording (e.g., 60s recording: 0.1 = 6 s)
pk_init_range = [5 40];                 % pk_init_range: start range for peaks (s)
pk_dur_range = [5 40];                  % pk_dur_range: duration range for peaks (s)                                     
pk_cf_range = [3 35];                   % pk_cf_range: peak center frequency range (Hz)    
pk_amp_range = [0.6 1.6];               % pk_amp_range: peak amplitude range (a.u.)  
pk_sd_range = [1 2];                    % pk_sd_range: peak standard deviation range (Hz)  
min_pk_dist = 2.5;                      % min_pk_dist: minimum (frequency) distance between peaks (in largest sd of given peaks)
max_n_peaks = 4;                        % max_n_peaks: maximum number of peaks across simulation
n = 10;                                 % n: number of simulations
fs = 200;                               % fs: sampling rate

os = gen_rand_TF_par_rot(60,ap_init_ranges,ap_chng_range,[ap_chng_start_range;ap_chng_dur_range],...
    pk_init_range,pk_dur_range,pk_cf_range,pk_amp_range,pk_sd_range,min_pk_dist,max_n_peaks,n,fs);
xs2 = zeros(n,len*fs);

for s = 1:n

    xs2(s,:) = sim_combined_mat_full_rot(len,fs,os(s).ap_init,os(s).pk_cf,...
        os(s).pk_amp,os(s).pk_sd,os(s).pk_times,[0 0.3 0.8].*ones(size(os(s).pk_times,1),1),...
        [0 0.3 0.8].*ones(size(os(s).pk_times,1),1),...
        [os(s).ap_chng os(s).ap_dyn_prop],os(s).ap_rotf,1);
    % NOTE: [0 0.3 0.8]s are vestigial (left from previous dynamic peak components, see challenge 1)
end



%% Functions to generate simulated time series
function outStruct = gen_rand_TF_par_rot(len,ap_ranges,dyn_ap_ranges,ap_dyn_time_range,...
    pk_start_range, pk_len_range,pk_cf_range,pk_amp_range,pk_sd_range,min_pk_dist,max_n_peaks,...
    n,tune_fs)

% Generate simulation parameters for sim_combined_mat_full_rot from uniform
% distributions
% len                   : simulation length (in seconds)
% ap_ranges             : 2x2 [exp_min exp_max; rotf_min rotf_max]
% dyn_ap_ranges         : 1x2 [exp_chng_min exp_chng_max]
% ap_dyn_time_range     : 2x2 [ap_chng_start_prop_min ap_chng_start_prop_max; ap_chng_dur_prop_min ap_chng_dur_prop_max]
% pk_start_range        : 1x2 [min_start_time max_start_time]
% pk_len_range          : 1x2 [pk_len_min pk_len_max]
% pk_cf_range           : 1x2 [pk_cf_min pk_cf_max]
% pk_amp_range          : 1x2 [pk_amp_min pk_amp_max]
% pk_sd_range           : 1x2 [pk_sd_min pk_sd_max]
% min_pk_dist           : Minimum distance from other peaks, in sd of other peaks
% max_n_peaks           : Maximum number of unique peaks in a simulation
% n                     : Number of unique simulations
% tune_fs               : frequency to tune generated parameters to (in Hz)
% -- -- -- -- 
% outstruct             : output structure, containing simulation
%                         parameters for sim_combined_mat_full_rot

    outStruct = struct();

    for s = 1:n
        % construct aperiodic parameters
        ap_rotf = 10.^(log10(ap_ranges(2,1))+diff(log10(ap_ranges(2,:)).*rand));
        ap_init = tune_param([ap_ranges(1,1)+diff(ap_ranges(1,:).*rand) -5.9],tune_fs);
        ap_init(2) = -5.9+ap_init(1).*log10(ap_rotf);
        ap_chng = [tune_param(dyn_ap_ranges(1)+diff(dyn_ap_ranges.*rand),tune_fs) 0];
        ap_chng(2) = ap_chng(1).*log10(ap_rotf);
            ap_dyn_start = ap_dyn_time_range(1,1)+diff(ap_dyn_time_range(1,:).*rand);
        ap_dyn_prop = [ap_dyn_start min([0.95 ap_dyn_start+ap_dyn_time_range(2,1)+diff(ap_dyn_time_range(2,:)).*rand])];
        
        % Determine number of peaks
        npeaks = floor((max_n_peaks+1).*rand);
        if npeaks > max_n_peaks % The off-chance of pulling rand = 1
            npeaks = npeaks-1;
        end
        % Generate peak parameters
        pk_times = zeros(npeaks,2); pk_cf = zeros(npeaks,1); pk_amp = zeros(npeaks,1); pk_sd = zeros(npeaks,1); 
        int_timevec = zeros(npeaks,2);
        for pk = 1:npeaks
            pk_time_profile = [pk_start_range(1)+diff(pk_start_range).*rand pk_len_range(1)+diff(pk_len_range).*rand];
            ctr = 1;
            while sum(pk_time_profile) > len
                pk_time_profile = [pk_start_range(1)+diff(pk_start_range).*rand pk_len_range(1)+diff(pk_len_range).*rand];
                ctr = ctr + 1;
                if ctr > 100
                   error('Randomization process failed: peak time profile') 
                end
            end
            % Detect if peaks overlap temporally, then ensure that they 
            % do not overlap in frequency space.
            if pk == 1
                pk_times(pk,:) = pk_time_profile;
                pk_cf(pk) = pk_cf_range(1)+diff(pk_cf_range).*rand;
                pk_amp(pk) = pk_amp_range(1)+diff(pk_amp_range).*rand;
                pk_sd(pk) = pk_sd_range(1)+diff(pk_sd_range).*rand;
                int_timevec(pk,:) = [pk_time_profile(1) sum(pk_time_profile)];
            else
                pk_times(pk,:) = pk_time_profile;
                pk_amp(pk) = pk_amp_range(1)+diff(pk_amp_range).*rand;
                pk_sd(pk) = pk_sd_range(1)+diff(pk_sd_range).*rand;
                % determine whether(/which) peaks overlap temporally.
                overlaps_t = zeros(npeaks,1);
                grc = 0.5;
                for p = 1:pk-1
                    overlaps_t(p) = ((int_timevec(p,1)-grc) < pk_time_profile(1) && pk_time_profile(1) < (int_timevec(p,2))+grc) || ...
                        ((int_timevec(p,1)-grc) < sum(pk_time_profile) && sum(pk_time_profile) < (int_timevec(p,2))+grc) || ...
                        ((pk_time_profile(1)-grc) < int_timevec(p,1) && int_timevec(p,1) < (sum(pk_time_profile)+grc));
                end
                overlaps_t = logical(overlaps_t);
                pk_cf(pk) = pk_cf_range(1)+diff(pk_cf_range).*rand;
                if any(overlaps_t)
                    % generate a safe zone of frequencies
                    sd = max(pk_sd(overlaps_t),pk_sd(pk));
                    zn = [pk_cf(overlaps_t)-min_pk_dist.*sd pk_cf(overlaps_t)+min_pk_dist.*sd];
                    ctr = 1;
                    overlaps_f = (zn(:,1) < pk_cf(pk) & pk_cf(pk) < zn(:,2));
                    % If there are overlapping peaks, relocate them
                    while any(overlaps_f)
                        ctr = ctr+1;
                        if ctr > 100
                            error('Randomization process failed: peaks overlap') 
                        end
                        pk_cf(pk) = pk_cf_range(1)+diff(pk_cf_range).*rand;
                        overlaps_f = (zn(:,1) < pk_cf(pk) & pk_cf(pk) < zn(:,2));
                    end                   
                end
                int_timevec(pk,:) = [pk_time_profile(1) sum(pk_time_profile)];
            end
        end
        
        outStruct(s).ap_init = ap_init;
        outStruct(s).ap_chng = ap_chng;
        outStruct(s).ap_dyn_prop = ap_dyn_prop;
        outStruct(s).ap_rotf = ap_rotf;
        outStruct(s).pk_times = pk_times;
        outStruct(s).pk_cf = pk_cf;
        outStruct(s).pk_amp = pk_amp;
        outStruct(s).pk_sd = pk_sd;
        
    end
end

function param = tune_param(param, fs)
    % Rounds parameters to sampling rate
    param = ceil(param.*fs)./fs;
    
end

function xs = sim_combined_mat_full_rot(varargin)
% A comprehensive signal simulator 
% Based on works in NeuroDSP/sim by Cole et al. (2019)
% https://github.com/neurodsp-tools/neurodsp/blob/main/paper/paper.md
% doi:10.21105/joss.01272
% varargin{1} = len (s)
% varargin{2} = fs (Hz)
% varargin{3} = exponent (-2 for example)
% varargin{4} = [cf ...] in Hz
% varargin{5} = [hgt ...] in log10(Hz)
% varargin{6} = [bw ...] in Hz
% varargin{10} = [dsl, fracStart, fracStop]
% varargin{11} = maintain offset (1: yes, 0: no)
% varargin{12} = number of simulations

    time = varargin{1};
    fs = varargin{2};
    ap_pars = varargin{3};
    
    freqs = linspace(0,fs./2,round(time.*fs./2));
    times = linspace(1./fs,time.*fs./fs,time.*fs);

    if length(varargin) > 3 % if a peak was provided
        cf = varargin{4};
        hgt = varargin{5};
        stdev = varargin{6};
        tps = ceil(varargin{7}.*fs)./fs;
        if length(varargin) > 7
            cfr = varargin{8};
            sdr = varargin{9};
            ap_dyn = varargin{10};
            rot_f = varargin{11};
        end
        nS = varargin{12};
        
        % Ensures proportional representation of frequencies 
        % (done here to speed up simulations of the same signal)
        cos_norms = zeros(size(freqs));
        for f = 1:length(freqs)
            cos_norms(f) = norm(cos(2.*pi.*freqs(f).*times),2).^2;
        end

        rel_heights = zeros(length(times),length(freqs));
        for peak = 1:length(cf)
            envMat = zeros(1,time*fs)';
            % Generate centre frequency shift vector
            cfrMat = zeros(size(times))'; 
            if length(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,2)+1:tps(peak,1)*fs+ ...
                tps(peak,2)*fs*cfr(peak,3)) ~= round(round(cfr(peak,3)-cfr(peak,2),2)*tps(peak,2)*fs)
                cfrMat(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,2):tps(peak,1)*fs+ ...
                    tps(peak,2)*fs*cfr(peak,3))=linspace(0,cfr(peak,1),...
                    round(round(cfr(peak,3)-cfr(peak,2),2)*tps(peak,2)*fs));
            else
                cfrMat(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,2)+1:tps(peak,1)*fs+ ...
                    tps(peak,2)*fs*cfr(peak,3))=linspace(0,cfr(peak,1),...
                    round(round(cfr(peak,3)-cfr(peak,2),2)*tps(peak,2)*fs));
            end
            cfrMat(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,3):end) = cfr(peak,1);
            % Generate standard deviation shift vector
            sdrMat = zeros(size(times))';
            if length(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,2)+1:tps(peak,1)*fs+ ...
                tps(peak,2)*fs*sdr(peak,3)) ~= round(round(sdr(peak,3)-sdr(peak,2),2)*tps(peak,2)*fs)
                sdrMat(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,2):tps(peak,1)*fs+ ...
                    tps(peak,2)*fs*sdr(peak,3))=linspace(0,sdr(peak,1),round((sdr(peak,3)-sdr(peak,2))*tps(peak,2)*fs));
            else
                sdrMat(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,2)+1:tps(peak,1)*fs+ ...
                    tps(peak,2)*fs*sdr(peak,3))=linspace(0,sdr(peak,1),round((sdr(peak,3)-sdr(peak,2))*tps(peak,2)*fs));
            end
            sdrMat(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,3):end) = sdr(peak,1);
            % Generate amplitude vector (with Tukey kernel)
            if length(tps(peak,1)*fs+1:(tps(peak,1)+tps(peak,2))*fs) ~= round(tps(peak,2)*fs)
                envMat(tps(peak,1)*fs:(tps(peak,1)+tps(peak,2))*fs) = tukeywin(tps(peak,2)*fs,0.4);
            else
                envMat(tps(peak,1)*fs+1:(tps(peak,1)+tps(peak,2))*fs) = tukeywin(tps(peak,2).*fs,0.4);
            end
            % Generate peak power profile across time
            rh_tmp = ones(size(freqs)).*(envMat*hgt(peak)).*...
                exp(-(ones(size(envMat)).*freqs-ones(size(freqs)).*...
                (cf(peak)+cfrMat)).^2./(2.*(ones(size(freqs)).*(stdev(peak)+sdrMat).^2)));
            % Combine with previous peaks
            rel_heights = rel_heights +rh_tmp;
        end
        % Generate white noise time series
        xs_tmp = sim_aperiodic_mat_mult(time*fs,0,1,nS);
        % Generate simulated neural time series
        xs = sim_spec_full_rot(xs_tmp, ap_pars, fs, ap_dyn, rot_f, cos_norms, rel_heights);
    else
        xs = sim_aperiodic_mat_mult(time*fs,slope,1,1);
    end
    
end

function xs = sim_aperiodic_mat_mult(points,exp,powNoise,n)

    % Simulate colored noise
    % points    : number of samples
    % exp       : aperiodic exponent (range: [-2 2])
    % powNoise  : can be used to change offset from -5.9 (in white noise)
    % nS        : number of simulations
    % -- -- -- --
    % xs        : simulated aperiodic time series
    
noise = dsp.ColoredNoise(exp,points,n);
ns = noise()';

xs = (ns*powNoise);
end

function xs = sim_spec_full_rot(aper_xs, ap_pars, fs, ap_dyn, rot_f, cos_norms, rel_heights)
    
    % Simulate neural time series
    % Based on works in NeuroDSP/sim by Cole et al. (2019)
    % https://github.com/neurodsp-tools/neurodsp/blob/main/paper/paper.md
    % doi:10.21105/joss.01272
    % aper_xs       : white noise time series (same #samples as desired time series)
    %                   NOTE: as is, offset initialized at -5.9 a.u.
    % ap_pars       : initial aperiodic parameters (only position 1 is necessary)
    % fs            : sampling rate (Hz)
    % ap_dyn        : [change_in_exponent change_in_offset start_of_change(in prop) end_of_change(in prop)]
    % rot_f         : rotation frequency (Hz; determines offset, proportional to exponent)
    % cos_norms     : ensures proportional representation of each frequency
    % rel_heights   : periodic "spectrogram"
    % -- -- -- --
    % xs            : simulated neural time series

    slen = length(aper_xs);
    times = linspace(1./fs,slen./fs,slen);
    
    % generate aperiodic fft
    sig_fft = fft(aper_xs);
    sig_fft = sig_fft(1:round(end./2))./100;
    freqs = linspace(0,fs./2,round(length(aper_xs)./2));
    % Calculate exponent at all time samples
    slrMat = ones(size(times))'.*ap_pars(1); % before aperiodic shift
    % during aperiodic shift
    if length(ap_dyn(3)*slen+1:ap_dyn(4)*slen) ~= round(ap_dyn(4).*slen-(ap_dyn(3).*slen))
        slrMat(ap_dyn(3)*slen:ap_dyn(4)*slen)=ap_pars(1)+linspace(0,ap_dyn(1),...
            round(ap_dyn(4).*slen-(ap_dyn(3).*slen)));
    else
        slrMat(ap_dyn(3)*slen+1:ap_dyn(4)*slen)=ap_pars(1)+linspace(0,ap_dyn(1),...
            round(ap_dyn(4).*slen-(ap_dyn(3).*slen)));
    end
    % after aperiodic shift
    slrMat(ap_dyn(4)*slen+1:end) = ap_dyn(1)+ap_pars(1);
    
    % Generate aperiodic "spectrogram"
    apMat = ones(size(rel_heights)); % initialize
    apMat(:,2:end) = ones(length(times),1)*sig_fft(2:end).*(freqs(2:end)./rot_f).^(-slrMat./2);    
    apDiff = -real(apMat);
    apDiff(apDiff > 0) = 0;
    % Combine aperiodic and periodic "spectrograms"
    cos_coeffs = (-real(apMat) +apDiff +sqrt(real(apMat).^2 + ...
        (10.^rel_heights - 1).*abs(apMat).^2))./cos_norms; 
    % Generate simulated neural time series
    xs = sum(cos_coeffs.*cos(times'*2*pi*freqs+2*pi*rand(size(freqs))),2,'omitnan')';
end

function xs = sim_combined_mat5(varargin)
% A comprehensive signal simulator
% Based on works in NeuroDSP/sim by Cole et al. (2019)
% https://github.com/neurodsp-tools/neurodsp/blob/main/paper/paper.md
% doi:10.21105/joss.01272
% varargin{1} = time (s)
% varargin{2} = fs (Hz)
% varargin{3} = slope (-2 for example)
% varargin{4} = [cf ...] in Hz
% varargin{5} = [hgt ...] in log10(Hz)
% varargin{6} = [bw ...] in Hz
% varargin{10} = [dsl, fracStart, fracStop]
% varargin{11} = maintain offset (1: yes, 0: no)
% varargin{12} = number of simulations

    time = varargin{1};
    fs = varargin{2};
    slope = varargin{3};
    
    freqs = linspace(0,fs./2,round(time.*fs./2));
    times = linspace(1./fs,time.*fs./fs,time.*fs);
    
    

    if length(varargin) > 3 % if a peak was provided
        cf = varargin{4};
        hgt = varargin{5};
        stdev = varargin{6};
        tps = varargin{7};

        if length(varargin) > 7
            cfr = varargin{8};
            sdr = varargin{9};
            slr = varargin{10};
            slr = [slr(1) 0 slr(2:3)];
            rot_f = varargin{11};
        end
        
        
        n = varargin{12};
        xs = zeros(n,time*fs);
        
        cos_norms = zeros(size(freqs));
        for f = 1:length(freqs)
            cos_norms(f) = norm(cos(2.*pi.*freqs(f).*times),2).^2;
        end
        tic
        rel_heights = zeros(length(times),length(freqs));
        for peak = 1:length(cf)
            envMat = zeros(1,time*fs)';
            % Generate centre frequency shift vector
            cfrMat = zeros(size(times))'; 
            cfrMat(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,2)+1:tps(peak,1)*fs+ ...
                tps(peak,2)*fs*cfr(peak,3))=linspace(0,cfr(peak,1),...
                round(cfr(peak,3)-cfr(peak,2),2)*tps(peak,2)*fs);
            cfrMat(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,3):end) = cfr(peak,1);
            % Generate standard deviation shift vector
            sdrMat = zeros(size(times))';
            sdrMat(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,2)+1:tps(peak,1)*fs+ ...
                tps(peak,2)*fs*sdr(peak,3))=linspace(0,sdr(peak,1),(sdr(peak,3)-sdr(peak,2))*tps(peak,2)*fs);
            sdrMat(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,3):end) = sdr(peak,1);
            % Generate amplitude envelope
            envMat(tps(peak,1)*fs+1:(tps(peak,1)+tps(peak,2))*fs) = tukeywin(tps(peak,2).*fs,0.4);
            % Generate peak "spectrogram"
            rh_tmp = ones(size(freqs)).*(envMat*hgt(peak)).*...
                exp(-(ones(size(envMat)).*freqs-ones(size(freqs)).*...
                (cf(peak)+cfrMat)).^2./(2.*(ones(size(freqs)).*(stdev(peak)+sdrMat).^2)));
            % Combine with other peaks
            rel_heights = rel_heights +rh_tmp;
        end
        % Generate white noise spectrum
        xst = sim_aperiodic_mat_mult(time*fs,0,1,n);
        for s = 1:n
            % Combine aperiodic and periodic components
            xs(s,:) = sim_spec_full_rot(xst(s,:), slope, fs, slr, rot_f, cos_norms, rel_heights);
            if ~mod(s,n./10)
                disp(['trial ' num2str(s) ' of ' num2str(n)])
            end
        end
    else
        xs = sim_aperiodic_mat_mult(time*fs,slope,1,1);
    end

end

function outStruct = gen_rand_TF_par_rot_kneerev(ap_ranges,dyn_ap_ranges,ap_dyn_time_range,...
    pk_start_range,pk_end_range,pk_amp_range,pk_sd_range,...
    n,tune_fs)

% Generate simulation parameters for sim_combined_mat_full_rot from uniform
% distributions
% len                   : simulation length (in seconds)
% ap_ranges             : 2x2 [exp_min exp_max; rotf_min rotf_max]
% dyn_ap_ranges         : 1x2 [exp_chng_min exp_chng_max]
% ap_dyn_time_range     : 2x2 [ap_chng_start_prop_min ap_chng_start_prop_max; ap_chng_dur_prop_min ap_chng_dur_prop_max]
% pk_start_range        : 1x2 [min_start_time max_start_time]
% pk_len_range          : 1x2 [pk_len_min pk_len_max]
% pk_cf_range           : 1x2 [pk_cf_min pk_cf_max]
% pk_amp_range          : 1x2 [pk_amp_min pk_amp_max]
% pk_sd_range           : 1x2 [pk_sd_min pk_sd_max]
% min_pk_dist           : Minimum distance from other peaks, in sd of other peaks
% max_n_peaks           : Maximum number of unique peaks in a simulation
% n                     : Number of unique simulations
% tune_fs               : frequency to tune generated parameters to (in Hz)
% -- -- -- -- 
% outstruct             : output structure, containing simulation
%                         parameters for sim_combined_mat_full_rot

    outStruct = struct();

    for s = 1:n
        % construct aperiodic parameters
        ap_rotf = 10.^(log10(ap_ranges(2,1))+diff(log10(ap_ranges(2,:)).*rand));
        ap_init = tune_param([ap_ranges(1,1)+diff(ap_ranges(1,:).*rand) -5.9],tune_fs);
        ap_init(2) = -5.9+ap_init(1).*log10(ap_rotf);
        ap_knee = 30.*rand; % Add knee parameter
        ap_chng = [tune_param(dyn_ap_ranges(1)+diff(dyn_ap_ranges.*rand),tune_fs) 0];
        ap_chng(2) = ap_chng(1).*log10(ap_rotf);
            ap_dyn_start = ap_dyn_time_range(1,1)+diff(ap_dyn_time_range(1,:).*rand);
        ap_dyn_prop = [ap_dyn_start min([0.95 ap_dyn_start+ap_dyn_time_range(2,1)+diff(ap_dyn_time_range(2,:)).*rand])];
        
        % Determine number of peaks
        npeaks = 2;
        % Generate peak parameters
        pk_times = zeros(npeaks,2); pk_cf = zeros(npeaks,1); pk_amp = zeros(npeaks,1); pk_sd = zeros(npeaks,1); 
        for pk = 1:npeaks
            pk_time_profile = [pk_start_range(1)+diff(pk_start_range).*rand pk_end_range(1)+diff(pk_end_range).*rand];
            pk_time_profile(2) = diff(pk_time_profile);
            if pk == 1
                pk_times(pk,:) = pk_time_profile;
                pk_cf(pk) = 3+27*rand;
                pk_amp(pk) = pk_amp_range(1)+diff(pk_amp_range).*rand;
                pk_sd(pk) = pk_sd_range(1)+diff(pk_sd_range).*rand;
            else
                pk_times(pk,:) = pk_time_profile;
                pk_amp(pk) = pk_amp_range(1)+diff(pk_amp_range).*rand;
                pk_sd(pk) = pk_sd_range(1)+diff(pk_sd_range).*rand;
                % determine whether(/which) peaks overlap temporally.
                pk_cf(pk) = max([pk_cf(1)+max(pk_sd).*2.5 30])+(80-max([pk_cf(1)+max(pk_sd).*2.5 30])).*rand;
            end
        end
        
        outStruct(s).ap_init = ap_init;
        outStruct(s).ap_chng = ap_chng;
        outStruct(s).ap_dyn_prop = ap_dyn_prop;
        outStruct(s).ap_rotf = ap_rotf;
        outStruct(s).ap_knee = ap_knee;
        outStruct(s).pk_times = pk_times;
        outStruct(s).pk_cf = pk_cf;
        outStruct(s).pk_amp = pk_amp;
        outStruct(s).pk_sd = pk_sd;
        
    end
end

function xs = sim_combined_mat_full_rot_knee(varargin)
% A comprehensive signal simulator, but this time with an aperiodic knee!
% Based on works in NeuroDSP/sim by Cole et al. (2019)
% https://github.com/neurodsp-tools/neurodsp/blob/main/paper/paper.md
% doi:10.21105/joss.01272
% varargin{1} = len (s)
% varargin{2} = fs (Hz)
% varargin{3} = exponent (-2 for example)
% varargin{4} = [cf ...] in Hz
% varargin{5} = [hgt ...] in log10(Hz)
% varargin{6} = [bw ...] in Hz
% varargin{7} = [start dur ...] for each simulated peak (in s)
% varargin{8} = [cfr ...] in Hz
% varargin{9} = [sdr ...] in Hz
% varargin{10} = [dsl, fracStart, fracStop]
% varargin{11} = aperiodic rotation frequency in Hz
% varargin{12} = aperiodic knee frequency
% varargin{13} = number of simulations (1, loop outside fn)

    time = varargin{1};
    fs = varargin{2};
    ap_pars = varargin{3};
    
    freqs = linspace(0,fs./2,round(time.*fs./2));
    times = linspace(1./fs,time.*fs./fs,time.*fs);

    if length(varargin) > 3 % if a peak was provided
        cf = varargin{4};
        hgt = varargin{5};
        stdev = varargin{6};
        tps = ceil(varargin{7}.*fs)./fs;
        if length(varargin) > 7
            cfr = varargin{8};
            sdr = varargin{9};
            ap_dyn = varargin{10};
            rot_f = varargin{11};
            knee = varargin{12};
        end
        nS = varargin{13};
        
        % Ensures proportional representation of frequencies 
        % (done here to speed up simulations of the same signal)
        cos_norms = zeros(size(freqs));
        for f = 1:length(freqs)
            cos_norms(f) = norm(cos(2.*pi.*freqs(f).*times),2).^2;
        end

        rel_heights = zeros(length(times),length(freqs));
        for peak = 1:length(cf)
            envMat = zeros(1,time*fs)';
            % Generate centre frequency shift vector
            cfrMat = zeros(size(times))'; 
            if length(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,2)+1:tps(peak,1)*fs+ ...
                tps(peak,2)*fs*cfr(peak,3)) ~= round(round(cfr(peak,3)-cfr(peak,2),2)*tps(peak,2)*fs)
                cfrMat(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,2):tps(peak,1)*fs+ ...
                    tps(peak,2)*fs*cfr(peak,3))=linspace(0,cfr(peak,1),...
                    round(round(cfr(peak,3)-cfr(peak,2),2)*tps(peak,2)*fs));
            else
                cfrMat(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,2)+1:tps(peak,1)*fs+ ...
                    tps(peak,2)*fs*cfr(peak,3))=linspace(0,cfr(peak,1),...
                    round(round(cfr(peak,3)-cfr(peak,2),2)*tps(peak,2)*fs));
            end
            cfrMat(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,3):end) = cfr(peak,1);
            % Generate standard deviation shift vector
            sdrMat = zeros(size(times))';
            if length(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,2)+1:tps(peak,1)*fs+ ...
                tps(peak,2)*fs*sdr(peak,3)) ~= round(round(sdr(peak,3)-sdr(peak,2),2)*tps(peak,2)*fs)
                sdrMat(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,2):tps(peak,1)*fs+ ...
                    tps(peak,2)*fs*sdr(peak,3))=linspace(0,sdr(peak,1),round((sdr(peak,3)-sdr(peak,2))*tps(peak,2)*fs));
            else
                sdrMat(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,2)+1:tps(peak,1)*fs+ ...
                    tps(peak,2)*fs*sdr(peak,3))=linspace(0,sdr(peak,1),round((sdr(peak,3)-sdr(peak,2))*tps(peak,2)*fs));
            end
            sdrMat(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,3):end) = sdr(peak,1);
            % Generate amplitude vector (with Tukey kernel)
            if length(tps(peak,1)*fs+1:(tps(peak,1)+tps(peak,2))*fs) ~= round(tps(peak,2)*fs)
                envMat(tps(peak,1)*fs:(tps(peak,1)+tps(peak,2))*fs) = tukeywin(tps(peak,2)*fs,0.4);
            else
                envMat(tps(peak,1)*fs+1:(tps(peak,1)+tps(peak,2))*fs) = tukeywin(tps(peak,2).*fs,0.4);
            end
            % Generate peak power profile across time
            rh_tmp = ones(size(freqs)).*(envMat*hgt(peak)).*...
                exp(-(ones(size(envMat)).*freqs-ones(size(freqs)).*...
                (cf(peak)+cfrMat)).^2./(2.*(ones(size(freqs)).*(stdev(peak)+sdrMat).^2)));
            % Combine with previous peaks
            rel_heights = rel_heights +rh_tmp;
        end
        % Generate white noise time series
        xs_tmp = sim_aperiodic_mat_mult(time*fs,0,1,nS);
        % Generate simulated neural time series with aperiodic knee
        xs = sim_spec_full_rot_knee(xs_tmp, ap_pars, fs, ap_dyn, rot_f, knee, cos_norms, rel_heights);
    else
        xs = sim_aperiodic_mat_mult(time*fs,slope,1,1);
    end
    
end

function xs = sim_spec_full_rot_knee(aper_xs, ap_pars, fs, ap_dyn, rot_f, knee_f, cos_norms, rel_heights)
    
    % Simulate neural time series, this time with an aperiodic knee!
    % Based on works in NeuroDSP/sim by Cole et al. (2019)
    % https://github.com/neurodsp-tools/neurodsp/blob/main/paper/paper.md
    % doi:10.21105/joss.01272
    % aper_xs       : white noise time series (same #samples as desired time series)
    %                   NOTE: as is, offset initialized at -5.9 a.u.
    % ap_pars       : initial aperiodic parameters (only position 1 is necessary)
    % fs            : sampling rate (Hz)
    % ap_dyn        : [change_in_exponent change_in_offset start_of_change(in prop) end_of_change(in prop)]
    % rot_f         : rotation frequency (Hz; determines offset, proportional to exponent)
    % knee_f        : aperiodic knee frequency (Hz)
    % cos_norms     : ensures proportional representation of each frequency
    % rel_heights   : periodic "spectrogram"
    % -- -- -- --
    % xs            : simulated neural time series

    slen = length(aper_xs);
    times = linspace(1./fs,slen./fs,slen);
    
    % generate aperiodic fft
    sig_fft = fft(aper_xs);
    sig_fft = sig_fft(1:round(end./2))./100;
    freqs = linspace(0,fs./2,round(length(aper_xs)./2));
    % Calculate exponent at all time samples
    slrMat = ones(size(times))'.*ap_pars(1); % before aperiodic shift
    % during aperiodic shift
    if length(ap_dyn(3)*slen+1:ap_dyn(4)*slen) ~= round(ap_dyn(4).*slen-(ap_dyn(3).*slen))
        slrMat(ap_dyn(3)*slen:ap_dyn(4)*slen)=ap_pars(1)+linspace(0,ap_dyn(1),...
            round(ap_dyn(4).*slen-(ap_dyn(3).*slen)));
    else
        slrMat(ap_dyn(3)*slen+1:ap_dyn(4)*slen)=ap_pars(1)+linspace(0,ap_dyn(1),...
            round(ap_dyn(4).*slen-(ap_dyn(3).*slen)));
    end
    % after aperiodic shift
    slrMat(ap_dyn(4)*slen+1:end) = ap_dyn(1)+ap_pars(1);
    
    % Generate aperiodic "spectrogram"
    apMat = ones(size(rel_heights)); % initialize
    apMat(:,2:end) = ones(length(times),1)*sig_fft(2:end).*(freqs(2:end)./rot_f).^(-slrMat./2);
    
    % Knee compensation
    apMat(:,2:end) = apMat(:,2:end)./sqrt((knee_f./freqs(2:end).^slrMat+1));
    
    apDiff = -real(apMat);
    apDiff(apDiff > 0) = 0;
    % Combine aperiodic and periodic "spectrograms"
    cos_coeffs = (-real(apMat) +apDiff +sqrt(real(apMat).^2 + ...
        (10.^rel_heights - 1).*abs(apMat).^2))./cos_norms; 
    % Generate simulated neural time series
    xs = sum(cos_coeffs.*cos(times'*2*pi*freqs+2*pi*rand(size(freqs))),2,'omitnan')';
end
