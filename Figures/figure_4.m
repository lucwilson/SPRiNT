%% Figure 4 - LEMON SPRiNT data

% Load data
maxdims = 0;
mindims = 1000;
work = whos('SPRiNT*'); % files are named SPRiNT001, SPRiNT002, etc.
% Initialize table columns
subject = cell(numel(work),1);
eo_ap_slope_mean = nan(numel(work),1);
eo_ap_slope_std = nan(numel(work),1);
eo_ap_off_mean = nan(numel(work),1);
eo_ap_off_std = nan(numel(work),1);
eo_alpha_cf_mean = nan(numel(work),1);
eo_alpha_cf_std = nan(numel(work),1);
eo_alpha_pow_mean = nan(numel(work),1);
eo_alpha_pow_std = nan(numel(work),1);

ec_ap_slope_mean = nan(numel(work),1);
ec_ap_slope_std = nan(numel(work),1);
ec_ap_off_mean = nan(numel(work),1);
ec_ap_off_std = nan(numel(work),1);
ec_alpha_cf_mean = nan(numel(work),1);
ec_alpha_cf_std = nan(numel(work),1);
ec_alpha_pow_mean = nan(numel(work),1);
ec_alpha_pow_std = nan(numel(work),1);

eo_R2 = nan(numel(work),1);
ec_R2 = nan(numel(work),1);
glob_R2 = nan(numel(work),1);

eo_MAE = nan(numel(work),1);
ec_MAE = nan(numel(work),1);
glob_MAE = nan(numel(work),1);

% Determine dimensions
for t = 1:numel(work)
    baseline = eval(work(t).name);
    maxdims = max([size(squeeze(baseline.SPRiNT.peak_models(1,:,:)),1) maxdims]);
    mindims = min([size(squeeze(baseline.SPRiNT.peak_models(1,:,:)),1) mindims]);
end

% Initialize spectrogram matrices
stft_spec = zeros(mindims,40);
SPRiNT_spec = zeros(mindims,40);
ts = baseline.Time(1:mindims);

for t = 1:numel(work)
    % Two naming schemes were used
    baseline = eval(work(t).name);
    if strcmp('sub',baseline.DataFile(1:3))
        subject{t} = baseline.DataFile(1:10);
    else
        subject{t} = baseline.DataFile(6:15);
    end
    for chan = 1:length(baseline.SPRiNT.channel)
        % Identify channel of interest
        if strcmp(baseline.SPRiNT.channel(chan).name,'Oz')
            eo = [119:238,359:478]; % Eyes open times
            ec = [1:118,239:358,479:min([598,length(baseline.Time)])]; % Eyes closed times
            
            % Extract exponent descriptives
            ap_exp = [baseline.SPRiNT.channel(chan).aperiodics.exponent];
            eo_ap_slope_mean(t) = mean(ap_exp(eo));
            ec_ap_slope_mean(t) = mean(ap_exp(ec));
            eo_ap_slope_std(t) = std(ap_exp(eo));
            ec_ap_slope_std(t) = std(ap_exp(ec));
            
            % Extract offset descriptives
            ap_off = [baseline.SPRiNT.channel(chan).aperiodics.offset];
            eo_ap_off_mean(t) = mean(ap_off(eo));
            ec_ap_off_mean(t) = mean(ap_off(ec));
            eo_ap_off_std(t) = std(ap_off(eo));
            ec_ap_off_std(t) = std(ap_off(ec));
            
            % Extract alpha peak descriptives
            alp_pk = baseline.SPRiNT.channel(chan).peaks([baseline.SPRiNT.channel(chan).peaks.center_frequency]<=14 & [baseline.SPRiNT.channel(chan).peaks.center_frequency]>=6);
            [~, eo_tmp, ~] = intersect([alp_pk.time],baseline.Time(eo));
            [~, ec_tmp, ~] = intersect([alp_pk.time],baseline.Time(ec));
            eo_alpha_cf_mean(t) = mean([alp_pk(eo_tmp).center_frequency]);
            ec_alpha_cf_mean(t) = mean([alp_pk(ec_tmp).center_frequency]);
            eo_alpha_cf_std(t) = std([alp_pk(eo_tmp).center_frequency]);
            ec_alpha_cf_std(t) = std([alp_pk(ec_tmp).center_frequency]);
            eo_alpha_pow_mean(t) = mean([alp_pk(eo_tmp).amplitude]);
            eo_alpha_pow_std(t) = std([alp_pk(eo_tmp).amplitude]);
            ec_alpha_pow_mean(t) = mean([alp_pk(ec_tmp).amplitude]);
            ec_alpha_pow_std(t) = std([alp_pk(ec_tmp).amplitude]);
            
            % extract goodness of fit stats
            r2_eo = corrcoef(log10(baseline.TF(chan,eo,:)),log10(baseline.SPRiNT.SPRiNT_models(chan,eo,:)));
            eo_R2(t) = r2_eo(2).^2;
            eo_MAE(t) = mean(mean(abs(log10(baseline.TF(chan,eo,:))-log10(baseline.SPRiNT.SPRiNT_models(chan,eo,:)))));
            
            r2_ec = corrcoef(log10(baseline.TF(chan,ec,:)),log10(baseline.SPRiNT.SPRiNT_models(chan,ec,:)));
            ec_R2(t) = r2_ec(2).^2;
            ec_MAE(t) = mean(mean(abs(log10(baseline.TF(chan,ec,:))-log10(baseline.SPRiNT.SPRiNT_models(chan,ec,:)))));
            
            r2_full = corrcoef(log10(baseline.TF(chan,:,:)),log10(baseline.SPRiNT.SPRiNT_models(chan,:,:)));
            glob_R2(t) = r2_full(2).^2;
            glob_MAE(t) = mean(mean(abs(log10(baseline.TF(chan,:,:))-log10(baseline.SPRiNT.SPRiNT_models(chan,:,:)))));
            % Generate mean spectra for Figure 4
            stft_spec = stft_spec + log10(squeeze(baseline.TF(chan,1:mindims,:)));
            SPRiNT_spec = SPRiNT_spec + log10(squeeze(baseline.SPRiNT.SPRiNT_models(chan,1:mindims,:)));
            break
        else
            continue
        end
    end
    
    
end

% Organize data into table
Oz = table(subject,eo_ap_slope_mean,eo_ap_slope_std,eo_ap_off_mean,eo_ap_off_std,eo_alpha_cf_mean,eo_alpha_cf_std,eo_alpha_pow_mean,eo_alpha_pow_std,...
    ec_ap_slope_mean,ec_ap_slope_std,ec_ap_off_mean,ec_ap_off_std,ec_alpha_cf_mean,ec_alpha_cf_std,ec_alpha_pow_mean,ec_alpha_pow_std,eo_R2,eo_MAE,ec_R2,ec_MAE,glob_R2,glob_MAE);
Oz = sortrows(Oz,1);
% Write table, if desired
% writetable(Oz,'SPRiNT_Oz.csv','Delimiter',','); 

% Average spectrograms by number of samples
stft_spec = stft_spec./t;
SPRiNT_spec = SPRiNT_spec./t;

% Generate Figure 4
figure('Position',[500 500 1200 800])

% Panel A: psd+specparam models for eyes-closed (ec) and eyes-open (eo) resting state.
subplot(3,2,1), hold on
plot(1:40,specparam_ec,'--b')
plot(1:40,psd_ec,'b')
plot(1:40,specparam_eo,'--y')
plot(1:40,psd_eo,'y')
xlim([1 40])
ylim([-13.5 -11])
xlabel('Frequency (Hz)','FontSize',12)
ylabel('Log power (a.u.)','FontSize',12)

% Panel C: stft for sub016, with model spectra adjacent
subplot(3,2,3), hold on
imagesc(1.5:0.5:size(sub016_spec,1)/2+1,1:40,sub016_spec');
ax = gca;
ax.YDir = 'normal';
ylim([1 40])
xlim([1.5 size(SPRiNT_spec,1)/2+1])
c = colorbar;
c.Label.String = 'Log Power (a.u.)';
caxis([-14.1 -10.9])
xticks(60:60:300)
xlabel('Time (s)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)

subplot(3,2,4), hold on
plot(1:40,specparam016,'--b')
plot(1:40,psd016,'b')
plot(1:40,SPRiNT016,'--k')
plot(1:40,stft016,'k')
xlim([1 40])
ylim([-14.1 -10.9])
xlabel('Frequency (Hz)','FontSize',12)
ylabel('Log power (a.u.)','FontSize',12)

% Panel D: average SPRiNT model spectrogram, with s067 alpha cf
% distribution adjacent
subplot(3,2,5), hold on
imagesc(1.5:0.5:size(SPRiNT_spec,1)/2+1,1:40,SPRiNT_spec');
ax = gca;
ax.YDir = 'normal';
ylim([1 40])
xlim([1.5 size(SPRiNT_spec,1)/2+1])
c = colorbar;
c.Label.String = 'Log Power (a.u.)';
caxis([-14.1 -10.9])
xticks(60:60:300)
xlabel('Time (s)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)

subplot(3,2,6), hold on 
histogram(sub067_alpha,'BinEdges',6:0.5:14) % but axes flipped
