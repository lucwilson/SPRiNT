%% Data for figure 3

% Load data
work = whos('SPRiNT*'); % Data is saved as SPRiNT01, SPRiNT02, etc.
if numel(work) > 1
    sd = eval(work(1).name).SPRiNT;
    for t = 2:numel(work)
        sdtmp = eval(work(t).name).SPRiNT;
        sd.channel(end+1:end+length(sdtmp.channel)) = sdtmp.channel;
    end
else
    sd = eval(work(1).name).SPRiNT;
end

%
ts = [sd.channel(1).data.time];
len = 60; % length of recording, in seconds
fs = 200;
tvec = 1/fs:1/fs:60;
tuk = tukeywin(1000,0.4);

sdtmp = sd; % SPRiNT data structures
ostmp = os; % Simulation parameters
vstmp = struct(); % Expected output statistics (under perfect conditions)
rstmp = struct(); % SPRiNT output statistics

% Generate expected peak statistics
for sim = 1:length(sdtmp.channel)
    vstmp(sim).ap_exp = zeros(fs.*len,1);
    vstmp(sim).ap_off = zeros(fs.*len,1);
    
    vstmp(sim).ap_exp(1:round(ostmp(sim).ap_dyn_prop(1).*len.*fs)) = ostmp(sim).ap_init(1);
    vstmp(sim).ap_exp(round(ostmp(sim).ap_dyn_prop(1).*len.*fs)+1:round(ostmp(sim).ap_dyn_prop(2).*len.*fs)) = ...
        linspace(ostmp(sim).ap_init(1),ostmp(sim).ap_init(1)+ostmp(sim).ap_chng(1),round(ostmp(sim).ap_dyn_prop(2).*len.*fs)-round(ostmp(sim).ap_dyn_prop(1).*len.*fs));
    vstmp(sim).ap_exp(round(ostmp(sim).ap_dyn_prop(2).*len.*fs)+1:end) = ostmp(sim).ap_init(1)+ostmp(sim).ap_chng(1);
    vstmp(sim).ap_exp = spline(tvec,vstmp(sim).ap_exp,ts);
    
    vstmp(sim).ap_off(1:round(ostmp(sim).ap_dyn_prop(1).*len.*fs)) = ostmp(sim).ap_init(2);
    vstmp(sim).ap_off(round(ostmp(sim).ap_dyn_prop(1).*len.*fs)+1:round(ostmp(sim).ap_dyn_prop(2).*len.*fs)) = ...
        linspace(ostmp(sim).ap_init(2),ostmp(sim).ap_init(2)+ostmp(sim).ap_chng(2),round(ostmp(sim).ap_dyn_prop(2).*len.*fs)-round(ostmp(sim).ap_dyn_prop(1).*len.*fs));
    vstmp(sim).ap_off(round(ostmp(sim).ap_dyn_prop(2).*len.*fs)+1:end) = ostmp(sim).ap_init(2)+ostmp(sim).ap_chng(2);
    vstmp(sim).ap_off = spline(tvec,vstmp(sim).ap_off,ts);
    
    vstmp(sim).pk_cf = nan(size(ostmp(sim).pk_cf,1),length(ts));
    vstmp(sim).pk_amp = nan(size(ostmp(sim).pk_cf,1),length(ts));
    vstmp(sim).pk_sd = nan(size(ostmp(sim).pk_cf,1),length(ts));
    vstmp(sim).pk_exist = zeros(size(ostmp(sim).pk_cf,1),length(ts));
    
    for pk = 1:size(ostmp(sim).pk_cf,1)
        vstmp(sim).pk_exist(pk,:) = (ostmp(sim).pk_times(pk,1) <= ts) & (ts <= ostmp(sim).pk_times(pk,1)+ostmp(sim).pk_times(pk,2));
        vstmp(sim).pk_cf(pk,logical(vstmp(sim).pk_exist(pk,:))) = ostmp(sim).pk_cf(pk);
        vstmp(sim).pk_amp(pk,logical(vstmp(sim).pk_exist(pk,:))) = spline(1:1000,tuk,linspace(1,1000,sum(vstmp(sim).pk_exist(pk,:)))).*ostmp(sim).pk_amp(pk);
        vstmp(sim).pk_sd(pk,logical(vstmp(sim).pk_exist(pk,:))) = ostmp(sim).pk_sd(pk);
    end
    
    rstmp(sim).n_peaks = size(ostmp(sim).pk_cf,1);
    
end

% Extract peak statistics
sdtmp = sd;
for sim = 1:length(sdtmp.channel)
    npeak_vec = zeros(1,length(ts));
    npeak_exp = sum(~isnan(vstmp(sim).pk_cf),1);
    if isempty(npeak_exp)
        npeak_exp = zeros(1,length(ts));
    end
    peak_props = zeros(5,7,6);
    for i = 1:length(ts)
        npeak_vec(i) = sum(round([sdtmp.channel(sim).peaks.time],1)==round(ts(i),1));
    end
    for i = 1:length(ts)
        pks_tmp = [sdtmp.channel(sim).peaks(round([sdtmp.channel(sim).peaks.time],1)==round(ts(i),1)).center_frequency];
        for j = 1:length(pks_tmp)
            if pks_tmp(j) < 8
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,1) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,1)+1;
            elseif 8 <= pks_tmp(j) && pks_tmp(j) < 13
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,2) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,2)+1;
            elseif 13 <= pks_tmp(j) && pks_tmp(j) < 18
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,3) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,3)+1;
            elseif 18 <= pks_tmp(j) && pks_tmp(j) < 23
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,4) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,4)+1;
            elseif 23 <= pks_tmp(j) && pks_tmp(j) < 28
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,5) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,5)+1;
            elseif 28 <= pks_tmp(j) && pks_tmp(j) < 33
                peak_props(npeak_exp(i)+1,npeak_vec(i)+1,6) = peak_props(npeak_exp(i)+1,npeak_vec(i)+1,6)+1;
            else
            end
        end
    end
    peaks = sdtmp.channel(sim).peaks;
    keep = zeros(1,length(peaks));
    zn = [ostmp(sim).pk_cf-2.*ostmp(sim).pk_sd ostmp(sim).pk_cf+2.*ostmp(sim).pk_sd];
    if size(ostmp(sim).pk_cf,1)
        for pk = 1:size(ostmp(sim).pk_cf,1)
            keep = keep | ((ostmp(sim).pk_times(pk,1) <= [peaks.time]) & ([peaks.time] <= ostmp(sim).pk_times(pk,1)+ostmp(sim).pk_times(pk,2)) & (zn(pk,1) <= [peaks.center_frequency]) & (zn(pk,2) >= [peaks.center_frequency]));
        end
    else
        keep = zeros(1,length(peaks));
    end
    discard = ~keep;
    rstmp(sim).resid_peaks = peaks(discard);
    sdtmp.channel(sim).peaks(discard) = [];
    peaks(discard) = [];    
    if size(ostmp(sim).pk_cf,1)
        [tmp, order] = sort(ostmp(sim).pk_amp,'descend');
        rstmp(sim).exp_peaks = struct();
        for ind = 1:length(order)
            pk = order(ind);
            expi = find((ostmp(sim).pk_times(pk,1) <= [peaks.time]) &...
                ([peaks.time] <=...
                ostmp(sim).pk_times(pk,1)+ostmp(sim).pk_times(pk,2)) &...
                (zn(pk,1) <= [peaks.center_frequency]) & (zn(pk,2) >=...
                [peaks.center_frequency])); % filter by cf and time range
%             expi = find((zn(pk,1) <= [peaks.center_frequency]) & (zn(pk,2) >= [peaks.center_frequency])); % filter by cf range only
            exp_peaks = peaks(expi);
            if length(exp_peaks) > length(unique([exp_peaks.time]))
                for i = unique([exp_peaks.time])
                    pks = find([exp_peaks.time] == i);
                    if length(pks) >1 
                        [tmp, rmv] = max([exp_peaks([exp_peaks.time] == i).amplitude]);
                        pks(rmv) = [];
                        exp_peaks(pks) = [];
                        expi(pks) = [];
                    end   
                end
            end
            rstmp(sim).exp_peaks(pk).peak = peaks(expi);
            peaks(expi) = [];
        end
        if ~isempty(peaks)
            rstmp(sim).resid_peaks(length(rstmp(sim).resid_peaks)+1:length(rstmp(sim).resid_peaks)+length(peaks)) = peaks;
        end
    end

    rstmp(sim).ap_exp = [sdtmp.channel(sim).aperiodics.exponent];
    rstmp(sim).ap_off = [sdtmp.channel(sim).aperiodics.offset];
    rstmp(sim).pk_cf = nan(size(vstmp(sim).pk_cf));
    rstmp(sim).pk_amp = nan(size(vstmp(sim).pk_amp));
    rstmp(sim).pk_sd = nan(size(vstmp(sim).pk_sd));
    rstmp(sim).pk_det = nan(size(vstmp(sim).pk_cf));
    
    for pk = 1:size(ostmp(sim).pk_cf,1)
        [~, iMat, iPk] = intersect(ts,[rstmp(sim).exp_peaks(pk).peak.time]);
        rstmp(sim).pk_cf(pk,iMat) = [rstmp(sim).exp_peaks(pk).peak(iPk).center_frequency];
        rstmp(sim).pk_amp(pk,iMat) = [rstmp(sim).exp_peaks(pk).peak(iPk).amplitude];
        rstmp(sim).pk_sd(pk,iMat) = [rstmp(sim).exp_peaks(pk).peak(iPk).st_dev];
        rstmp(sim).pk_det(pk,~isnan(vstmp(sim).pk_cf(pk,:))) = 0;
        rstmp(sim).pk_det(pk,iMat) = 1;
    end
    
    rstmp(sim).expmod_peaks = [npeak_exp;npeak_vec];
    rstmp(sim).peak_props = peak_props;
    
end
%%
% Calculate sensitivity and specificity
sens_mat = zeros(length(rstmp),2);
spec_mat = zeros(length(rstmp),1);
for i = 1:length(rstmp)
    sens_mat(i,1) = sum(vstmp(i).pk_exist(:));
    sens_mat(i,2) = sum(rstmp(i).pk_det(logical(vstmp(i).pk_exist)));
    spec_mat(i) = length(rstmp(i).resid_peaks);
end
sens = sum(sens_mat(:,2))./sum(sens_mat(:,1));
spec = 1-sum(spec_mat)./length([sdtmp.channel(1).data.time])./length(sdtmp.channel);

% Extract a variety of peak descriptives and errors
err_timefreq = nan(40000,10);
c = 0;
for s = 1:length(sd.channel)
    for i = 1:size(vstmp(s).pk_cf,1)
        c = c+1;
        err_timefreq(c,1) = ostmp(s).pk_times(i,2); % duration, in s
        err_timefreq(c,2) = ostmp(s).pk_cf(i); % centre frequency
        err_timefreq(c,3) = mean(abs(rstmp(s).pk_cf(i,:)-vstmp(s).pk_cf(i,:)),'omitnan'); % cf MAE
        err_timefreq(c,4) = mean(abs(rstmp(s).pk_amp(i,:)-vstmp(s).pk_amp(i,:)),'omitnan'); % amp MAE
        err_timefreq(c,5) = mean(abs(rstmp(s).pk_sd(i,:)-vstmp(s).pk_sd(i,:)),'omitnan'); % sd MAE
        err_timefreq(c,6) = ostmp(s).pk_sd(i); % sd
        err_timefreq(c,7) = ostmp(s).pk_times(i,2).*ostmp(s).pk_cf(i); % ncycles
        err_timefreq(c,8) = sum(~isnan(rstmp(s).pk_cf(i,~isnan(vstmp(s).pk_cf(i,:)))))./sum(~isnan(vstmp(s).pk_cf(i,:))); % peak detection rate
        err_timefreq(c,9) = mean(vstmp(s).ap_exp(~isnan(vstmp(s).pk_cf(i,:)))); % mean aperiodic slope
        err_timefreq(c,10) = ostmp(s).pk_amp(i); % amplitude
    end
end
err_timefreq(c+1:end,:) = [];

% Calculate error by number of peaks
err_bypk = zeros(3,5);
err_bypk(1,:) = 0:4;
err_bypkmat = zeros(1150000,5);
c = 0;  
for s = 1:length(sd.channel)
    for i = 1:115
        c = c+1;
        if isempty(vstmp(s).pk_cf)
            err_bypkmat(c,1) = 0;
        else
            err_bypkmat(c,1) = sum(~isnan(vstmp(s).pk_cf(:,i)));
        end
        err_bypkmat(c,2) = sqrt(40.*sdtmp.channel(s).stats(i).MSE)./40; % MAE
    end
end
err_bypkmat(c+1:end,:) =[];
err_bypkcell = {err_bypkmat(err_bypkmat(:,1) == 0,2),err_bypkmat(err_bypkmat(:,1) == 1,2),...
    err_bypkmat(err_bypkmat(:,1) == 2,2),err_bypkmat(err_bypkmat(:,1) == 3,2),err_bypkmat(err_bypkmat(:,1) == 4,2)};

% calculate error for each peak parameter
err_peaks = zeros(400000,3);
err_aper = zeros(length(sd.channel).*length(sd.channel(1)),2);
ctr = 1; 
for s = 1:length(sd.channel)
    for i = 1:size(vstmp(s).pk_cf,1)
        a = abs(rstmp(s).pk_cf(i,:)-vstmp(s).pk_cf(i,:));
        b = abs(rstmp(s).pk_amp(i,:)-vstmp(s).pk_amp(i,:));
        c = abs(rstmp(s).pk_sd(i,:)-vstmp(s).pk_sd(i,:));
        idx = ~isnan(a);
        err_peaks(ctr:ctr+sum(idx)-1,1) = a(idx);
        err_peaks(ctr:ctr+sum(idx)-1,2) = b(idx);
        err_peaks(ctr:ctr+sum(idx)-1,3) = c(idx);
        ctr = ctr+sum(idx);
    end
    err_aper((s-1)*length(ts)+1:s*length(ts),1) = abs(rstmp(s).ap_off-vstmp(s).ap_off);
    err_aper((s-1)*length(ts)+1:s*length(ts),2) = abs(rstmp(s).ap_exp-vstmp(s).ap_exp);
end
err_peaks(ctr:end,:) = [];
err_peaks = log10(err_peaks);
err_aper = log10(err_aper);
mAE_cell = {err_aper(:,2),err_aper(:,1),err_peaks(:,1),err_peaks(:,2),err_peaks(:,3)}; % Exponent, offset, cf, amp,sd

npeak_mat = zeros(5,7);
for s = 1:length(sdtmp.channel)
    for i = 1:length([rstmp(s).expmod_peaks])
        npeak_mat(rstmp(s).expmod_peaks(1,i)+1,rstmp(s).expmod_peaks(2,i)+1) =...
            npeak_mat(rstmp(s).expmod_peaks(1,i)+1,rstmp(s).expmod_peaks(2,i)+1)+1;
    end
end
pp = zeros(5,7,6);
for i = 1:length(sdtmp.channel)
    pp = pp + rstmp(i).peak_props;
end
% Focus on first three freq ranges
pp(:,:,4:6) = [];

% Generate figure
figure('Position',[250 250 1000 900])
% Figure panel A: sample spectrogram
subplot(3,3,1), hold on
tf1 = pcolor(ts,tFOOOF01.tFOOOF.freqs,log10(squeeze(sd.tFOOOF_models(1,:,:)))');
tf1.EdgeColor = 'None';
xlim([ts(1) ts(end)])
ylim([1 40])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
caxis([-11.7 -7.2])
for p = 1:length(ostmp(1).pk_cf)
    pt = ostmp(1).pk_times(p,:);
    pcf = ostmp(1).pk_cf(p);
    psd = ostmp(1).pk_sd(p);
    plot([pt(1) pt(1) [pt(1) pt(1)]+pt(2) pt(1)],...
        [-psd psd psd -psd -psd]+pcf,'k')
    plot([sum(pt) ts(end)],[pcf pcf],'-.k')
end
plot([13.5 13.5 14.5 14.5 13.5],[1 40 40 1 1], 'r')
pos = get(gca, 'Position');
pos(3) = pos(3)-0.08;
pos(4) = pos(4)-0.12;
pos(2) = pos(2)+0.09+0.05;
set(gca, 'Position', pos)
text(pos(1)+59, pos(2)+29, ['Peak 1: [' num2str(round(os(1).pk_cf(2),1)) ', '...
    num2str(round(os(1).pk_amp(2),1)) ', ' num2str(round(os(1).pk_sd(2),1)) ']'],'FontSize',8)
text(pos(1)+59, pos(2)+12, ['Peak 2: [' num2str(round(os(1).pk_cf(1),1)) ', '...
    num2str(round(os(1).pk_amp(1),1)) ', ' num2str(round(os(1).pk_sd(1),1)) ']'],'FontSize',8)
xl = [pos(1) pos(1)+pos(3).*os(1).ap_dyn_prop(1)];
yl = [pos(2)+pos(4)+0.002 pos(2)+pos(4)+0.002];
annotation('line',xl,yl)
text(xl(1)+1,44,['[' num2str(round(os(1).ap_init(2),1)) ', '...
    num2str(round(os(1).ap_init(1),1)) ']'],'FontSize',8)
xl = [pos(1)+pos(3).*os(1).ap_dyn_prop(1) pos(1)+pos(3).*sum(os(1).ap_dyn_prop)];
yl = [pos(2)+pos(4)+0.002 pos(2)+pos(4)+0.012];
annotation('line',xl,yl)
xl = [pos(1)+pos(3).*sum(os(1).ap_dyn_prop) pos(1)+pos(3)];
yl = [pos(2)+pos(4)+0.012 pos(2)+pos(4)+0.012];
annotation('line',xl,yl)
text(pos(1)+44.5,48,['[' num2str(round(os(1).ap_init(2)+os(1).ap_chng(2),1)) ', '...
    num2str(round(os(1).ap_init(1)+os(1).ap_chng(1),1)) ']'],'FontSize',8)
text(pos(1)+24,50,'Aperiodic','FontSize',8)
xl = [0.1592 0.1592];
yl = [0.78 0.72]+0.05;
annotation('line',xl,yl,'Color','r')
xl = [0.1592 0.20];
yl = [0.72 0.72]+0.05;
annotation('arrow',xl,yl,'Color','r')
% This code replaces the spectrogram (above) with a time-resolved spectrum
% (below)
% % % plot(1:40,-vstmp(1).ap_exp(18).*log10(1:40)+vstmp(1).ap_off(18))
% plot(1:40,log10(squeeze(tFOOOF01.TF(1,26,:))),'k')
% plot(1:40,log10(tFOOOF01.tFOOOF.channel(1).data(26).fooofed_spectrum),'r')
% plot(1:40,-tFOOOF01.tFOOOF.channel(1).aperiodics(26).exponent.*log10(1:40)+tFOOOF01.tFOOOF.channel(1).aperiodics(26).offset,'--b')
% % legend('Original spectrum','Full model fit','Aperiodic fit','FontSize',8)
% xlabel('Frequency (Hz)')
% ylabel('Power')
% yticks([])
% xticks(0:10:40)
% pos = get(gca, 'Position');
% % pos(3) = pos(3)-0.12;
% % pos(1) = pos(1)+0.12;
% pos(3) = pos(3)-0.09;
% pos(1) = pos(1)+0.09;
% pos(4) = pos(4)-0.15+0.05;
% % pos(2) = pos(2);
% % pos(2) = pos(2)-0.15;
% set(gca, 'Position', pos)

% Panel B: error by spectral parameter
subplot(3,3,2)
swatch = [52 72 148]./255;
violin(mAE_cell,'facecolor',swatch,'facealpha',1);
yticks([-4 -3 -2 -1 0 1])
yticklabels({'10^{-4}','10^{-3}','10^{-2}','0.1','1','10'})
ylim([-3.2 1])
ylabel('|error|')
xticks(1:6)
ax = gca(); 
row1 = {'Aperiodic exponent','Aperiodic offset','Center frequency','Amplitude','Standard deviation'};
row2 = {'(a.u. Hz^{-1})','(a.u.)','(Hz)','(a.u.)','(Hz)'};
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
xticklabels(tickLabels); 
pos = get(gca, 'Position');
pos(4) = pos(4)+0.05;
pos(3) = pos(3)+0.282;
set(gca, 'Position', pos)

% Panel C1: Peak detection rate by peak center frequency
subplot(3,3,4), hold on
pars = [0.85 4];
q = 0;
swatch = [108 156 167;236 189 78;198 97 43;162 55 27;91 61 28;86 34 98]./255;
scatter(err_timefreq(err_timefreq(:,2)<8 & err_timefreq(:,1)>q,2),err_timefreq(err_timefreq(:,2)<8 & err_timefreq(:,1)>q,8),6,swatch(1,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
scatter(err_timefreq(abs(err_timefreq(:,2)-10.5)<2.5 & err_timefreq(:,1)>q,2),err_timefreq(abs(err_timefreq(:,2)-10.5)<2.5 & err_timefreq(:,1)>q,8),6,swatch(2,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
scatter(err_timefreq(abs(err_timefreq(:,2)-15.5)<2.5 & err_timefreq(:,1)>q,2),err_timefreq(abs(err_timefreq(:,2)-15.5)<2.5 & err_timefreq(:,1)>q,8),6,swatch(3,:),'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
scatter(err_timefreq(err_timefreq(:,2)>18 & err_timefreq(:,1)>q,2),err_timefreq(err_timefreq(:,2)>18 & err_timefreq(:,1)>q,8),6,[0.5 0.5 0.5],'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3)
xlabel('Center frequency (Hz)')
ylabel('p(peak)')
xlim([3 35])
pos = get(gca, 'Position');
pos(4) = pos(4)+0.03;
pos(2) = pos(2)-0.01;
set(gca, 'Position', pos)

% Panel C2: Peak detection rate by peak amplitude
subplot(3,3,7)
swatch = [1 26 121]./255;
scatter(err_timefreq(:,10),err_timefreq(:,8),6,swatch,'filled','MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
xlabel('Amplitude (a.u.)')
ylabel('p(peak)')
pos = get(gca, 'Position');
pos(4) = pos(4)+0.03;
pos(2) = pos(2)-0.02;
set(gca, 'Position', pos)

% Panel D: Number of recovered peaks by number of simulated peaks
subplot(3,3,8), hold on
violin(err_bypkcell,'facecolor',[52 72 148]./255,'facealpha',1);
yticks(0.02:0.02:0.08)
yticklabels(2:2:8)
ylim([0.012 0.088])
ylabel('Mean absolute error (x10^{-2})')
xlim([0.65 5.5])
xticks(1:5)
xticklabels(0:4);
xlabel('Number of simulated peaks')
pos = get(gca, 'Position');
pos(3) = pos(3)+0.282;
pos(4) = pos(4)-0.07;
pos(2) = pos(2)-0.02;
set(gca, 'Position', pos)

% Panel E: Error by number of peaks
subplot(3,3,5), hold on
swatch = [108 156 167;236 189 78;198 97 43;162 55 27;91 61 28;86 34 98]./255;
xs = [0 0 0 0 0 0 0 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 4];
ys = [0 1 2 3 4 5 6 0 1 2 3 4 5 6 0 1 2 3 4 5 6 0 1 2 3 4 5 6 0 1 2 3 4 5 6];
sz = 0.9*[npeak_mat(1,:)./sum(npeak_mat(1,:)) npeak_mat(2,:)./sum(npeak_mat(2,:)) npeak_mat(3,:)./sum(npeak_mat(3,:)) npeak_mat(4,:)./sum(npeak_mat(4,:)) npeak_mat(5,:)./sum(npeak_mat(5,:))];
xs = xs(logical(sz));
ys = ys(logical(sz));
sz = sz(logical(sz));
pos = get(gca, 'Position');
pos(3) = pos(3)+0.282;
pos(4) = pos(4)+0.14;
pos(2) = pos(2)-0.12;
rto = ((4.45./4.85)./(pos(3)./pos(4))).*(62./58);
for p = 1:length(xs)
    if ys(p) == 0
        draw_circ(xs(p),ys(p),sz(p),[52 72 148]./255,rto);
    elseif ys(p) > 4.5
        continue
    else
        draw_pie(xs(p),ys(p),sz(p),pp(xs(p)+1,ys(p)+1,:),swatch,rto);
    end
end
yticks(0:6)
xticks(0:4)
xlim([-0.35 4.5])
ylim([-0.35 4.1])
ylabel('Number of recovered peaks')
xlabel('Number of simulated peaks')
set(gca, 'Position', pos)
for i = 1:size(npeak_mat,1)
    text(i-1, 3.7, num2str(sum(npeak_mat(i,:))),'HorizontalAlignment','center')
end

%% Support functions
function draw_circ(x,y,rad,color,rto)
    xs = x+cos(linspace(0,2*pi,180)).*rad.*rto;
    ys = y+sin(linspace(0,2*pi,180)).*rad;
    patch(xs, ys, color,'EdgeColor','white');
end

function draw_pie(x,y,rad,props,palette,rto)

    tot = sum(props);
    init = pi/2;
    for fr = 1:length(props)
        pts_per = round(200*props(fr)/tot);
        ths = linspace(init,init-props(fr)/tot*2*pi,pts_per);
        if isempty(ths)
            continue
        end
        xs = [x x+ones(1,pts_per).*cos(ths).*rad.*rto];
        ys = [y y+ones(1,pts_per).*sin(ths).*rad];
        init = ths(end);
        patch(xs, ys, palette(fr,:),'EdgeColor','white');
        
    end
    plot([x x], [y y+rad],'w')
    
    
end

function[h,L,MX,MED,bw]=violin(Y,varargin)
% This function is modified from the following source:
% Hoffmann H, 2015: violin.m - Simple violin plot using matlab default kernel
% density estimation. INRES (University of Bonn), Katzenburgweg 5, 53115 Germany.
% hhoffmann@uni-bonn.de
% INPUT
%
% Y:     Data to be plotted, being either
%        a) n x m matrix. A 'violin' is plotted for each column m, OR
%        b) 1 x m Cellarry with elements being numerical colums of nx1 length.
%
% varargin:
% xlabel:    xlabel. Set either [] or in the form {'txt1','txt2','txt3',...}
% facecolor: FaceColor. (default [1 0.5 0]); Specify abbrev. or m x 3 matrix (e.g. [1 0 0])
% edgecolor: LineColor. (default 'k'); Specify abbrev. (e.g. 'k' for black); set either [],'' or 'none' if the mean should not be plotted
% facealpha: Alpha value (transparency). default: 0.5
% mc:        Color of the bars indicating the mean. (default 'k'); set either [],'' or 'none' if the mean should not be plotted
% medc:      Color of the bars indicating the median. (default 'r'); set either [],'' or 'none' if the mean should not be plotted
% bw:        Kernel bandwidth. (default []); prescribe if wanted as follows:
%            a) if bw is a single number, bw will be applied to all
%            columns or cells
%            b) if bw is an array of 1xm or mx1, bw(i) will be applied to cell or column (i).
%            c) if bw is empty (default []), the optimal bandwidth for
%            gaussian kernel is used (see Matlab documentation for
%            ksdensity()
%
% OUTPUT
%
% h:     figure handle
% L:     Legend handle
% MX:    Means of groups
% MED:   Medians of groups
% bw:    bandwidth of kernel

%defaults:
%_____________________
xL=[];
fc=[205 205 205]./255;
lc='k';
mc=[];%'k';
medc='r';
b=[]; %bandwidth
plotlegend=0;
plotmean=0;
plotmedian=1;
x = [];
%_____________________
%convert single columns to cells:
if iscell(Y)==0
    Y = num2cell(Y,1);
end
%get additional input parameters (varargin)
if isempty(find(strcmp(varargin,'xlabel')))==0
    xL = varargin{find(strcmp(varargin,'xlabel'))+1};
end
if isempty(find(strcmp(varargin,'facecolor')))==0
    fc = varargin{find(strcmp(varargin,'facecolor'))+1};
end
if isempty(find(strcmp(varargin,'edgecolor')))==0
    lc = varargin{find(strcmp(varargin,'edgecolor'))+1};
end
if isempty(find(strcmp(varargin,'facealpha')))==0
    alp = varargin{find(strcmp(varargin,'facealpha'))+1};
end
if isempty(find(strcmp(varargin,'mc')))==0
    if isempty(varargin{find(strcmp(varargin,'mc'))+1})==0
        mc = varargin{find(strcmp(varargin,'mc'))+1};
        plotmean = 1;
    else
        plotmean = 0;
    end
end
if isempty(find(strcmp(varargin,'medc')))==0
    if isempty(varargin{find(strcmp(varargin,'medc'))+1})==0
        medc = varargin{find(strcmp(varargin,'medc'))+1};
        plotmedian = 1;
    else
        plotmedian = 0;
    end
end
if isempty(find(strcmp(varargin,'bw')))==0
    b = varargin{find(strcmp(varargin,'bw'))+1};
    if length(b)==1
        disp(['same bandwidth bw = ',num2str(b),' used for all cols'])
        b=repmat(b,size(Y,2),1);
    elseif length(b)~=size(Y,2)
        warning('length(b)~=size(Y,2)')
        error('please provide only one bandwidth or an array of b with same length as columns in the data set')
    end
end
if isempty(find(strcmp(varargin,'plotlegend')))==0
    plotlegend = varargin{find(strcmp(varargin,'plotlegend'))+1};
end
if isempty(find(strcmp(varargin,'x')))==0
    x = varargin{find(strcmp(varargin,'x'))+1};
end
%%
if size(fc,1)==1
    fc=repmat(fc,size(Y,2),1);
end
%% Calculate the kernel density
i=1;
fa = linspace(1,1./size(Y,2),size(Y,2));
for i=1:size(Y,2)
    
    if isempty(b)==0
        [f, u, bb]=ksdensity(Y{i},linspace(min(Y{i}),max(Y{i}),800),'bandwidth',b(i));
    elseif isempty(b)
        [f, u, bb]=ksdensity(Y{i},linspace(min(Y{i}),max(Y{i}),800));
    end
    
    f=f/max(f)*0.3; %normalize
    F(:,i)=f;
    U(:,i)=u;
    MED(:,i)=nanmedian(Y{i});
    MX(:,i)=nanmean(Y{i});
    bw(:,i)=bb;
    prc25(:,i) = prctile(Y{i},25);
    prc75(:,i) = prctile(Y{i},75);
    IQR(:,i) = diff([prc25(:,i) prc75(:,i)]);
    
end
%%
%-------------------------------------------------------------------------
% Put the figure automatically on a second monitor
% mp = get(0, 'MonitorPositions');
% set(gcf,'Color','w','Position',[mp(end,1)+50 mp(end,2)+50 800 600])
%-------------------------------------------------------------------------
%Check x-value options
if isempty(x)
    x = zeros(size(Y,2));
    setX = 0;
else
    setX = 1;
    if isempty(xL)==0
        disp('_________________________________________________________________')
        warning('Function is not designed for x-axis specification with string label')
        warning('when providing x, xlabel can be set later anyway')
        error('please provide either x or xlabel. not both.')
    end
end
%% Plot the violins
i=1;
for i=i:size(Y,2)
    if isempty(lc) == 1
        if setX == 0
            h(i)=fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',fa(i),'EdgeColor','none');
        else
            h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',fa(i),'EdgeColor','none');
        end
    else
        if setX == 0
            h(i)=fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',fa(i),'EdgeColor',lc);
        else
            h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],fc(i,:),'FaceAlpha',fa(i),'EdgeColor',lc);
        end
    end
    hold on
    if setX == 0
        if plotmean == 1
            p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i)) ],[MX(:,i) MX(:,i)],mc,'LineWidth',2);
        end
        if plotmedian == 1
            plot(ones(2,1).*mean([interp1(U(:,i),F(:,i)+i,MED(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))]),...
                [max([min(U(:,i)) prc25(:,i)-1.5.*IQR(:,i)]) min([max(U(:,i)) prc75(:,i)+1.5.*IQR(:,i)])],'Color',[0.2 0.2 0.2])
            plot(ones(2,1).*mean([interp1(U(:,i),F(:,i)+i,MED(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))]),...
                [prc25(:,i) prc75(:,i)],'Color',[0.2 0.2 0.2],'LineWidth',7)
            scatter(mean([interp1(U(:,i),F(:,i)+i,MED(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))]),MED(:,i),...
                abs(diff([interp1(U(:,i),F(:,i)+i,MED(:,i)), interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))])).*40,[1 1 1],'Filled');
            plot([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],'k')            
        end
    elseif setX == 1
        if plotmean == 1
            p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i))+x(i)-i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i))+x(i)-i],[MX(:,i) MX(:,i)],mc,'LineWidth',2);
        end
        if plotmedian == 1
            p(2)=plot([interp1(U(:,i),F(:,i)+i,MED(:,i))+x(i)-i, interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))+x(i)-i],[MED(:,i) MED(:,i)],medc,'LineWidth',2);
        end
    end
end
%% Add legend if requested
if plotlegend==1 & plotmean==1 | plotlegend==1 & plotmedian==1
    
    if plotmean==1 & plotmedian==1
        L=legend([p(1) p(2)],'Mean','Median');
    elseif plotmean==0 & plotmedian==1
        L=legend([p(2)],'Median');
    elseif plotmean==1 & plotmedian==0
        L=legend([p(1)],'Mean');
    end
    
    set(L,'box','off','FontSize',14)
else
    L=[];
end
%% Set axis
if setX == 0
    axis([0.5 size(Y,2)+0.5, min(U(:)) max(U(:))]);
elseif setX == 1
    axis([min(x)-0.05*range(x) max(x)+0.05*range(x), min(U(:)) max(U(:))]);
end
%% Set x-labels
xL2={''};
i=1;
for i=1:size(xL,2)
    xL2=[xL2,xL{i},{''}];
end
% set(gca,'TickLength',[0 0],'FontSize',12)
box on
if isempty(xL)==0
    set(gca,'XtickLabel',xL2)
end
%-------------------------------------------------------------------------
end %of function
