%% Figure 4 - LEMON SPRiNT data

% Load data
maxdims = 0;
mindims = 1000;
work = whos('tFOOOF*'); % files are named tFOOOF001, tFOOOF002, etc.
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
    maxdims = max([size(squeeze(baseline.tFOOOF.peak_models(1,:,:)),1) maxdims]);
    mindims = min([size(squeeze(baseline.tFOOOF.peak_models(1,:,:)),1) mindims]);
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
    for chan = 1:length(baseline.tFOOOF.channel)
        % Identify channel of interest
        if strcmp(baseline.tFOOOF.channel(chan).name,'Oz')
            eo = [119:238,359:478]; % Eyes open times
            ec = [1:118,239:358,479:min([598,length(baseline.Time)])]; % Eyes closed times
            
            % Extract exponent descriptives
            ap_exp = [baseline.tFOOOF.channel(chan).aperiodics.exponent];
            eo_ap_slope_mean(t) = mean(ap_exp(eo));
            ec_ap_slope_mean(t) = mean(ap_exp(ec));
            eo_ap_slope_std(t) = std(ap_exp(eo));
            ec_ap_slope_std(t) = std(ap_exp(ec));
            
            % Extract offset descriptives
            ap_off = [baseline.tFOOOF.channel(chan).aperiodics.offset];
            eo_ap_off_mean(t) = mean(ap_off(eo));
            ec_ap_off_mean(t) = mean(ap_off(ec));
            eo_ap_off_std(t) = std(ap_off(eo));
            ec_ap_off_std(t) = std(ap_off(ec));
            
            % Extract alpha peak descriptives
            alp_pk = baseline.tFOOOF.channel(chan).peaks([baseline.tFOOOF.channel(chan).peaks.center_frequency]<=14 & [baseline.tFOOOF.channel(chan).peaks.center_frequency]>=6);
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
            r2_eo = corrcoef(log10(baseline.TF(chan,eo,:)),log10(baseline.tFOOOF.tFOOOF_models(chan,eo,:)));
            eo_R2(t) = r2_eo(2).^2;
            eo_MAE(t) = mean(mean(abs(log10(baseline.TF(chan,eo,:))-log10(baseline.tFOOOF.tFOOOF_models(chan,eo,:)))));
            
            r2_ec = corrcoef(log10(baseline.TF(chan,ec,:)),log10(baseline.tFOOOF.tFOOOF_models(chan,ec,:)));
            ec_R2(t) = r2_ec(2).^2;
            ec_MAE(t) = mean(mean(abs(log10(baseline.TF(chan,ec,:))-log10(baseline.tFOOOF.tFOOOF_models(chan,ec,:)))));
            
            r2_full = corrcoef(log10(baseline.TF(chan,:,:)),log10(baseline.tFOOOF.tFOOOF_models(chan,:,:)));
            glob_R2(t) = r2_full(2).^2;
            glob_MAE(t) = mean(mean(abs(log10(baseline.TF(chan,:,:))-log10(baseline.tFOOOF.tFOOOF_models(chan,:,:)))));
            % Generate mean spectra for Figure 4
            stft_spec = stft_spec + log10(squeeze(baseline.TF(chan,1:mindims,:)));
            SPRiNT_spec = SPRiNT_spec + log10(squeeze(baseline.tFOOOF.tFOOOF_models(chan,1:mindims,:)));
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
% writetable(Oz,'tFOOOF_Oz.csv','Delimiter',','); 

% Average spectrograms by number of samples
stft_spec = stft_spec./t;
SPRiNT_spec = SPRiNT_spec./t;

% Generate Figure 4
figure('Position',[500 500 900 800])

% Panel A: stft, followed by SPRiNT
subplot(3,1,1), hold on
imagesc(1.5:0.5:size(stft_spec)/2+1,1:40,stft_spec');
ax = gca;
ax.YDir = 'normal';
ylim([1 40])
xlim([1.5 size(stft_spec,1)/2+1])
c = colorbar;
c.Label.String = 'Log Power (a.u.)';
caxis([-14.1 -10.9])
xticks(60:60:300)
    pos = get(gca, 'Position');
    pos(2) = pos(2)+0.0;
    set(gca, 'Position', pos)
xlabel('Time (s)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)

o = 0.0;
alp = 0.8;
subplot(3,1,2), hold on
draw_box([0 60], [-1 1]+o, [1 26 121],alp);
draw_box([61.5 121.5], [-1 1]+o, [239 161 7],alp);
draw_box([123 183], [-1 1]+o, [1 26 121],alp);
draw_box([184.5 244.5], [-1 1]+o, [239 161 7],alp);
draw_box([246 300], [-1 1]+o, [1 26 121],alp);
ylim([-1.5 1.5])
xlim([1.5 299])
    pos = get(gca, 'Position');
    pos(2) = pos(2)+0.05;
    set(gca, 'Position', pos)
axis off

subplot(3,1,3),hold on
b = imagesc(1.5:0.5:size(stft_spec)/2+1,1:40,SPRiNT_spec');
ax = gca;
ax.YDir = 'normal';
ylim([1 40])
xlim([1.5 size(stft_spec,1)/2+1])
c = colorbar;
c.Label.String = 'Log Power (a.u.)';
caxis([-14.1 -10.9])
xticks(60:60:300)
xlabel('Time (s)','FontSize',12)
ylabel('Frequency (Hz)','FontSize',12)