%% Code for generating Figure 2

% Initialize wavelet-specparam matrices
ap_exp2 = zeros(10000,10601);
ap_off2 = zeros(10000,10601);
pk1err_cfr2 = zeros(10000,7365);
pk1err_amp2 = zeros(10000,7365);
pk1err_std2 = zeros(10000,7365);
pk1_det2 = zeros(10000,10601);
pk2err_cfr2 = zeros(10000,2001);
pk2err_amp2 = zeros(10000,2001);
pk2err_std2 = zeros(10000,2001);
pk2_det2 = zeros(10000,10601);
fooof_model2 = zeros(40,10601);
err_ap2 = zeros(10000,10601);
err_pk2 = zeros(10000,10601);
err_fm2 = zeros(10000,10601);

% Load data from wavelet-specparam
for b = 1:10
    % cd to each folder containing data for wavelet-specparam in 10% blocks
    cd(['MW' num2str(b)]) % YOUR APPROPRIATE DIRECTORY HERE %
    load('MW_stats.mat')
    load('MW_peak1.mat')
    load('MW_peak2.mat')
    ap_exp2((b-1)*1000+1:b*1000,:) = ap_exp;
    ap_off2((b-1)*1000+1:b*1000,:) = ap_off;

    fooof_model2 = fooof_model2 + fooof_model;

    err_ap2((b-1)*1000+1:b*1000,:) = err_ap;
    err_pk2((b-1)*1000+1:b*1000,:) = err_pk;
    err_fm2((b-1)*1000+1:b*1000,:) = err_fm;

    pk1err_cfr2((b-1)*1000+1:b*1000,:) = pk1err_cfr;
    pk1err_amp2((b-1)*1000+1:b*1000,:) = pk1err_amp;
    pk1err_std2((b-1)*1000+1:b*1000,:) = pk1err_std;
    pk1_det2((b-1)*1000+1:b*1000,:) = pk1_det;
    pk2err_cfr2((b-1)*1000+1:b*1000,:) = pk2err_cfr;
    pk2err_amp2((b-1)*1000+1:b*1000,:) = pk2err_amp;
    pk2err_std2((b-1)*1000+1:b*1000,:) = pk2err_std;
    pk2_det2((b-1)*1000+1:b*1000,:) = pk2_det;
end
ap_exp = ap_exp2;
ap_off = ap_off2;
fooof_model = fooof_model2./10;
err_ap = err_ap2;
err_pk = err_pk2;
err_fm = err_fm2;

pk1err_cfr = pk1err_cfr2;
pk1err_amp = pk1err_amp2;
pk1err_std = pk1err_std2;
pk2err_cfr = pk2err_cfr2;
pk2err_amp = pk2err_amp2;
pk2err_std = pk2err_std2;

err_cfr_pk1 = log10(abs(pk1err_cfr2-pk1(1,pk1_exist)));
err_amp_pk1 = log10(abs(pk1err_amp2-pk1(2,pk1_exist)));
err_std_pk1 = log10(abs(pk1err_std2-pk1(3,pk1_exist)));
pk1_det = pk1_det2;
err_cfr_pk2 = log10(abs(pk2err_cfr2-pk2(1,pk2_exist)));
err_amp_pk2 = log10(abs(pk2err_amp2-pk2(2,pk2_exist)));
err_std_pk2 = log10(abs(pk2err_std2-pk2(3,pk2_exist)));
pk2_det = pk2_det2;

clear ap_exp2 ap_off2 fooof_model2 err_ap2 err_pk2 err_fm2 pk1err_cfr2 pk1err_amp2 
clear pk1err_std2 pk2err_cfr2 pk2err_amp2 pk2err_std2 pk1_det2 pk2_det2
%%
% Load data from SPRiNT
cd(['*/SPRiNT_data' num2str(b)]) 
load('SPRiNT_output_fig2.mat')

% Generate expected spectrogram
mat = gen_peak_mat(60, 200, [8 18 8 8], [1.2 0.9 1.2 1.2], [1.2 1.4 1.2 1.2], [8 32;15 10;41 5;47 5],...
    [0 0.3 0.8;-3 0.3 0.7;0 0.3 0.8;0 0.3 0.8], [0 0.3 0.8;0 0.3 0.8;0 0.3 0.8;0 0.3 0.8],...
    1:40, 1.5:0.005:58.5);
mat2 = gen_ap_mat1(1.5,[-2.56 -1.41], [0.5 24 36], 1:40, 1.5:0.005:58.5);

matr = mat+mat2;

% Generate Figure 2 (est time: 8 mins)
c1 = [1 50 150]./255;
c2 = [239 161 7]./255;
figure('Position',[300 300 900 700])

% Ground truth spectrogram
subplot(3,3,1), hold on
xlabel('Time (s)')
ylabel('Frequency (Hz)')
yticks(10:10:40)
xticks(10:10:50)
xlim([1.5 58.5])
ylim([1 40])
imagesc(1.5:0.005:58.5,1:40,matr)
caxis([-5.55 -1.4]);

% Wavelet-specparam model spectrogram
subplot(3,3,2), hold on
xlabel('Time (s)')
ylabel('Frequency (Hz)')
yticks(10:10:40)
xticks(10:10:50)
xlim([1.5 58.5])
ylim([1 40])
imagesc(1.5:0.5:58.5,1:40,fooof_model)
caxis([-5.55 -1.4]);

% SPRiNT model spectrogram
subplot(3,3,3), hold on
xlabel('Time (s)')
ylabel('Frequency (Hz)')
yticks(10:10:40)
xticks(10:10:50)
xlim([3.5 56.5])
ylim([1 40])
imagesc(3.5:0.005:56.5,1:40,fooof_mean)
caxis([-5.55 -1.4]);

% Wavelet-specparam aperiodic parameters
subplot(3,3,4), hold on
colororder([1 50 150; 239 161 7]./255)
yyaxis left
xlabel('Time (s)')
ylabel('Exponent (Hz^{-1})')
yticks(1:0.25:2.5)
xticks(10:10:50)
xlim([1.5 58.5])
ylim([1.2 2.3])
patch([-20 -20 -21],[-20 -21 -20],c1,'FaceAlpha',0.4)
patch([-20 -20 -21],[-20 -21 -20],c2,'FaceAlpha',0.4)
plot([1.5 24 36 56.5],[-2.56 -2.56 -1.41 -1.41],'--k')
draw_2575(times,ap_exp,c1)
plot([1.5 24 36 56.5],[1.5 1.5 2 2],'--k')
yyaxis right
ylabel('Offset (a.u.)')
draw_2575(times,ap_off,c2)
xlim([1.5 58.5])
ylim([-3.25 -0.72])
plot([1.5 24 36 56.5],[-2.56 -2.56 -1.41 -1.41],'--k')
yyaxis left
legend('Exponent','Offset','Ground truth','Location','Northwest')
pos = get(gca, 'Position');
    pos(3) = pos(3)+0.12;
    set(gca, 'Position', pos) 
    
% SPRiNT aperiodic parameters
subplot(3,3,6), hold on
colororder([1 50 150; 239 161 7]./255)
yyaxis left
xlabel('Time (s)')
ylabel('Exponent (Hz^{-1})')
yticks(1:0.25:2.5)
xticks(10:10:50)
xlim([1.5 58.5])
ylim([1.2 2.3])
x2 = [times_s fliplr(times_s)];
inBetween = [prctile(ap_exp_s,25), fliplr(prctile(ap_exp_s,75))];
f = fill(x2, inBetween, c1);
f.EdgeColor = 'None';
f.FaceAlpha = 0.2;
plot(times_s, median(ap_exp_s,'omitnan'),'-','Color',c1)
plot([1.5 24 36 56.5],[1.5 1.5 2 2],'--k')
yyaxis right
ylabel('Offset (a.u.)')
xlim([1.5 58.5])
ylim([-3.25 -0.72])
x2 = [times_s fliplr(times_s)];
inBetween = [prctile(ap_off_s,25), fliplr(prctile(ap_off_s,75))];
f = fill(x2, inBetween, c2);
f.EdgeColor = 'None';
f.FaceAlpha = 0.2;
plot(times_s, median(ap_off_s,'omitnan'),'-','Color',c2)
plot([1.5 24 36 56.5],[-2.56 -2.56 -1.41 -1.41],'--k')
pos = get(gca, 'Position');
    pos(1) = pos(1)-0.12;
    pos(3) = pos(3)+0.12;
    set(gca, 'Position', pos) 

% Periodic peak errors
subplot(3,3,7), hold on
err_cfr_pk1s = log10(abs(pk1err_cfr_s-pk1_s(1,pk1_exist_s)));
err_amp_pk1s = log10(abs(pk1err_amp_s-pk1_s(2,pk1_exist_s)));
err_std_pk1s = log10(abs(pk1err_std_s-pk1_s(3,pk1_exist_s)));
tmp = abs(pk1err_cfr-pk1(1,pk1_exist)); tmp = tmp(:); tmp(tmp == 0) = 0.0001;
err_cfr_pk1 = log10(tmp); err_cfr_pk1 = err_cfr_pk1(:); err_cfr_pk1 = err_cfr_pk1(~isnan(err_cfr_pk1));
err_amp_pk1 = log10(abs(pk1err_amp-pk1(2,pk1_exist))); err_amp_pk1 = err_amp_pk1(:); err_amp_pk1 = err_amp_pk1(~isnan(err_amp_pk1));
err_std_pk1 = log10(abs(pk1err_std-pk1(3,pk1_exist))); err_std_pk1 = err_std_pk1(:); err_std_pk1 = err_std_pk1(~isnan(err_std_pk1));
err_cfr_pk2s = log10(abs(pk2err_cfr_s-pk2_s(1,pk2_exist_s)));
err_amp_pk2s = log10(abs(pk2err_amp_s-pk2_s(2,pk2_exist_s)));
err_std_pk2s = log10(abs(pk2err_std_s-pk2_s(3,pk2_exist_s)));
err_cfr_pk2 = log10(abs(pk2err_cfr-pk2(1,pk2_exist))); err_cfr_pk2 = err_cfr_pk2(:); err_cfr_pk2 = err_cfr_pk2(~isnan(err_cfr_pk2));
err_amp_pk2 = log10(abs(pk2err_amp-pk2(2,pk2_exist))); err_amp_pk2 = err_amp_pk2(:); err_amp_pk2 = err_amp_pk2(~isnan(err_amp_pk2));
err_std_pk2 = log10(abs(pk2err_std-pk2(3,pk2_exist))); err_std_pk2 = err_std_pk2(:); err_std_pk2 = err_std_pk2(~isnan(err_std_pk2));
fill([-20 -20 -21],[-20 -21 -20],c1,'FaceAlpha',0.4)
fill([-20 -20 -21],[-20 -21 -20],c2,'FaceAlpha',0.4)
% Wavelet-specparam
violin_plot(err_cfr_pk1,c1,1,0.4,0.4,1,1);
violin_plot(err_cfr_pk2,c2,1,0.4,0.4,1,0);
violin_plot(err_amp_pk1,c1,2,0.4,0.4,1,1);
violin_plot(err_amp_pk2,c2,2,0.4,0.4,1,0);
violin_plot(err_std_pk1,c1,3,0.7,0.4,1,1);
violin_plot(err_std_pk2,c2,3,0.9,0.4,1,0);
% SPRiNT
violin_plot(err_cfr_pk1s(:),c1,6,0.4,0.4,1,1);
violin_plot(err_cfr_pk2s(:),c2,6,0.4,0.4,1,0);
violin_plot(err_amp_pk1s(:),c1,7,0.4,0.4,1,1);
violin_plot(err_amp_pk2s(:),c2,7,0.4,0.4,1,0);
violin_plot(err_std_pk1s(:),c1,8,0.4,0.4,1,1);
violin_plot(err_std_pk2s(:),c2,8,0.4,0.4,1,0);
xlim([0 9])
ylim([-3 0.5])
yticks([-4 -3 -2 -1 0 1])
yticklabels({'10^{-4}','10^{-3}','10^{-2}','0.1','1','10'})
ylim([-3.2 0.5])
ylabel('|error|')
xticks([1 2 3 6 7 8])
row1 = {'Center frequency','Amplitude','Standard deviation','Center frequency','Amplitude','Standard deviation'};
row2 = {'(Hz)','(a.u.)','(Hz)','(Hz)','(a.u.)','(Hz)'};
c2a = [249 217 156]./255;
c1a = [153 173 213]./255;
labelArray = [row1; row2]; 
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
xticklabels(tickLabels); 
pos = get(gca, 'Position');
    pos(3) = pos(3)+0.56;
    set(gca, 'Position', pos) 
a = fill([0.35 0.5 0.5 0.35]-0.1,[-2.8 -2.8 -2.2 -2.2], [0.95 0.95 0.95]);
a.EdgeColor = 'none';
a = fill([0.65 0.8 0.8 0.65]-0.1,[-2.8 -2.8 -2.2 -2.2], [0.95 0.95 0.95]);
a.EdgeColor = 'none';
% Sensitivity and specificity bars
plot([0.35 0.35]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
plot([0.50 0.50]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
plot([0.65 0.65]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
fill([0.2 0.788 0.788 0.2]-0.1, [-2.24 -2.24 -2.46 -2.46],c1a);
text(0.73,-2.35,'98','HorizontalAlignment','left','VerticalAlignment','middle')
fill([0.2 0.602 0.602 0.2]-0.1, [-2.76 -2.76 -2.54 -2.54],c2a);
text(0.73,-2.65,'67','HorizontalAlignment','left','VerticalAlignment','middle')
plot([0.2 0.2]-0.1, [-2.8 -2.2],'-','Color',[0.4 0.4 0.4])
plot([0.8 0.8]-0.1, [-2.8 -2.2],'-','Color',[0.4 0.4 0.4])
text(0.5-0.1,-2.05,'Sensitivity','HorizontalAlignment','center','VerticalAlignment','middle')
a = fill([3.35 3.5 3.5 3.35]-0.1,[-2.8 -2.8 -2.2 -2.2], [0.95 0.95 0.95]);
a.EdgeColor = 'none';
a = fill([3.65 3.8 3.8 3.65]-0.1,[-2.8 -2.8 -2.2 -2.2], [0.95 0.95 0.95]);
a.EdgeColor = 'none';
plot([3.35 3.35]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
plot([3.50 3.50]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
plot([3.65 3.65]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
fill([3.2 3.578 3.578 3.2]-0.1, [-2.24 -2.24 -2.46 -2.46],c1a);
text(3.73,-2.35,'63','HorizontalAlignment','left','VerticalAlignment','middle')
fill([3.2 3.794 3.794 3.2]-0.1, [-2.76 -2.76 -2.54 -2.54],c2a);
text(3.73,-2.65,'99','HorizontalAlignment','left','VerticalAlignment','middle')
plot([3.2 3.2]-0.1, [-2.8 -2.2],'-','Color',[0.4 0.4 0.4])
plot([3.8 3.8]-0.1, [-2.8 -2.2],'-','Color',[0.4 0.4 0.4])
text(3.5-0.1,-2.05,'Specificity','HorizontalAlignment','center','VerticalAlignment','middle')
a = fill([5.35 5.5 5.5 5.35]-0.1,[-2.8 -2.8 -2.2 -2.2], [0.95 0.95 0.95]);
a.EdgeColor = 'none';
a = fill([5.65 5.8 5.8 5.65]-0.1,[-2.8 -2.8 -2.2 -2.2], [0.95 0.95 0.95]);
a.EdgeColor = 'none';
plot([5.35 5.35]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
plot([5.50 5.50]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
plot([5.65 5.65]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
fill([5.2 5.794 5.794 5.2]-0.1, [-2.24 -2.24 -2.46 -2.46],c1a);
text(5.73,-2.35,'99','HorizontalAlignment','left','VerticalAlignment','middle')
fill([5.2 5.770 5.770 5.2]-0.1, [-2.76 -2.76 -2.54 -2.54],c2a);
text(5.73,-2.65,'95','HorizontalAlignment','left','VerticalAlignment','middle')
plot([5.2 5.2]-0.1, [-2.8 -2.2],'-','Color',[0.4 0.4 0.4])
plot([5.8 5.8]-0.1, [-2.8 -2.2],'-','Color',[0.4 0.4 0.4])
text(5.4,-2.05,'Sensitivity','HorizontalAlignment','center','VerticalAlignment','middle')
a = fill([8.35 8.5 8.5 8.35]-0.1,[-2.8 -2.8 -2.2 -2.2], [0.95 0.95 0.95]);
a.EdgeColor = 'none';
a = fill([8.65 8.8 8.8 8.65]-0.1,[-2.8 -2.8 -2.2 -2.2], [0.95 0.95 0.95]);
a.EdgeColor = 'none';
plot([8.35 8.35]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
plot([8.50 8.50]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
plot([8.65 8.65]-0.1, [-2.8 -2.2],'-','Color',[0.5 0.5 0.5])
fill([8.2 8.776 8.776 8.2]-0.1, [-2.24 -2.24 -2.46 -2.46],c1a);
text(8.73,-2.35,'96','HorizontalAlignment','left','VerticalAlignment','middle')
fill([8.2 8.794 8.794 8.2]-0.1, [-2.76 -2.76 -2.54 -2.54],c2a);
text(8.73,-2.65,'99','HorizontalAlignment','left','VerticalAlignment','middle')
plot([8.2 8.2]-0.1, [-2.8 -2.2],'-','Color',[0.4 0.4 0.4])
plot([8.8 8.8]-0.1, [-2.8 -2.2],'-','Color',[0.4 0.4 0.4])
text(8.5-0.1,-2.05,'Specificity','HorizontalAlignment','center','VerticalAlignment','middle')
legend('Alpha peak','Beta peak','Location','North')

%% Representative code to generate data structures for Figure 2

work = whos('SPRiNT*'); % Load data structures as SPRiNT01, SPRiNT02, etc.

fooof_mean = zeros(40,115);
err_aps = zeros(10000,115);
err_pks = zeros(10000,115);
err_fms = zeros(10000,115);

mat = gen_peak_mat(60, 200, [8 18 8 8], [1.2 0.9 1.2 1.2], [1.2 1.4 1.2 1.2], [8 32;15 10;41 5;47 5],...
    [0 0.3 0.8;-3 0.3 0.7;0 0.3 0.8;0 0.3 0.8], [0 0.3 0.8;0 0.3 0.8;0 0.3 0.8;0 0.3 0.8],...
    1:40, 1.5:0.5:58.5);
mat2 = gen_ap_mat1(1.5,[-2.56 -1.41], [0.5 24 36], 1:40, 1.5:0.5:58.5);
matr = mat+mat2;

if numel(work) > 1
    sd2 = eval(work(1).name).SPRiNT;
    for t = 2:numel(work)
        sdtmp = eval(work(t).name).SPRiNT;
        fooof_mean = fooof_mean + squeeze(mean(log10(sdtmp.SPRiNT_models),1))';
        fooof_mean = fooof_mean + squeeze(mean(log10(sdtmp.SPRiNT_models),1))';
        sd2.channel(end+1:end+length(sdtmp.channel)) = sdtmp.channel;
        for chan = 1:length(sdtmp.channel)
            for d = 1:115
                err_aps((t-1)*1000+chan,d) = mean(squeeze(log10(sdtmp.aperiodic_models(chan,d,:))) - mat2(:,d));
                err_pks((t-1)*1000+chan,d) = mean(squeeze(log10(sdtmp.peak_models(chan,d,:))) - mat(:,d));
                err_fms((t-1)*1000+chan,d) = mean(squeeze(log10(sdtmp.SPRiNT_models(chan,d,:))) - matr(:,d));
            end
        end
    end
else
    sd2 = eval(work(1).name).SPRiNT;
end

times = [sd2.channel(1).data.time];

sdtmp = sd2;

pk1 = nan(3,115);
pk2 = nan(3,115);

tukey = tukeywin(10000,0.4);
[~, pk1_exist,~] = intersect(times,[8:0.5:40,41:0.5:46,47:0.5:52]);

pk1(1,pk1_exist) = 8;
pk1(3,pk1_exist) = 1.2;
[~, pk1_exist,~] = intersect(times,8:0.5:40);
pk1(2,pk1_exist) = spline(1:10000,tukey,round((1:length(pk1_exist))./length(pk1_exist).*10000,0)).*1.2;
[~, pk1_exist,~] = intersect(times,41:0.5:46);
pk1(2,pk1_exist) = spline(1:10000,tukey,round((1:length(pk1_exist))./length(pk1_exist).*10000,0)).*1.2;
[~, pk1_exist,~] = intersect(times,47:0.5:52);
pk1(2,pk1_exist) = spline(1:10000,tukey,round((1:length(pk1_exist))./length(pk1_exist).*10000,0)).*1.2;
[~, pk1_exist,~] = intersect(times,[8:0.5:40,41:0.5:46,47:0.5:52]);

[~, pk2_exist,~] = intersect(times,15:0.5:25);
pk2(3,pk2_exist) = 1.4;
pk2(2,pk2_exist) = spline(1:10000,tukey,round((1:length(pk2_exist))./length(pk2_exist).*10000,0)).*0.9;
pk2(1,pk2_exist) = 18;
[~, pk2_exist,~] = intersect(times,18:0.5:22);
pk2(1,pk2_exist) = pk2(1,pk2_exist)-linspace(0,3,length(pk2_exist));
[~, pk2_exist,~] = intersect(times,22:0.5:25);
pk2(1,pk2_exist) = 15;
[~, pk2_exist,~] = intersect(times,15:0.5:25);

pk1err_cfr = nan(length(sdtmp.channel),sum(~isnan(pk1(1,:))));
pk1err_amp = nan(length(sdtmp.channel),sum(~isnan(pk1(1,:))));
pk1err_std = nan(length(sdtmp.channel),sum(~isnan(pk1(1,:))));
pk1_det = zeros(length(sdtmp.channel),length(times));

pk2err_cfr = nan(length(sdtmp.channel),sum(~isnan(pk2(1,:))));
pk2err_amp = nan(length(sdtmp.channel),sum(~isnan(pk2(1,:))));
pk2err_std = nan(length(sdtmp.channel),sum(~isnan(pk2(1,:))));
pk2_det = zeros(length(sdtmp.channel),length(times));

pit = struct('discard',[]);
i = 1;
pk1_seen = 0;
pk2_seen = 0;

sdtmp1 = sdtmp;
sdtmp2 = sdtmp;

for chan = 1:length(sdtmp.channel)
    peaks = sdtmp.channel(chan).peaks;
    pks = logical([peaks.center_frequency]>20.5) | logical(abs([peaks.center_frequency]-12)<1.5) |...
        logical([peaks.center_frequency]< 5.5);
    pit(i).discard = sdtmp.channel(chan).peaks(pks);
    i = i+1;
    peaks = sdtmp.channel(chan).peaks(~pks);
    sdtmp1.channel(chan).peaks = peaks([peaks.center_frequency] < 11);
    sdtmp2.channel(chan).peaks = peaks([peaks.center_frequency] > 11);
    sdtmp.channel(chan).peaks(pks) = [];
end

for chan = 1:length(sdtmp.channel)
peaks1 = sdtmp1.channel(chan).peaks;
for pk = 1:length(peaks1)
    peaks1(pk).time = round(peaks1(pk).time,3);
end
if length(unique([peaks1.time])) ~= length([peaks1.time])
    for i = round(unique([peaks1.time]),3)
        pks = find([peaks1.time] == i);
        if length(pks) >1 
            [tmp, rmv] = max([peaks1([peaks1.time] == i).amplitude]);
            pks(rmv) = [];
            peaks1(pks) = [];
        end
    end
end
sdtmp1.channel(chan).peaks = peaks1;

peaks2 = sdtmp2.channel(chan).peaks;
for pk = 1:length(peaks2)
    peaks2(pk).time = round(peaks2(pk).time,3);
end
if length(unique([peaks2.time])) ~= length([peaks2.time])
    for i = round(unique([peaks2.time]),3)
        pks = find([peaks2.time] == i);
        if length(pks) >1 
            [tmp, rmv] = max([peaks2([peaks2.time] == i).amplitude]);
            pks(rmv) = [];
            peaks2(pks) = [];
        end
    end
end
sdtmp2.channel(chan).peaks = peaks2;

end

for chan = 1:length(sdtmp.channel)
    sdtmp.channel(chan).clustered_peaks = [];
    sdtmp.channel(chan).clustered_peaks(1).peaks = sdtmp1.channel(chan).peaks;
    sdtmp.channel(chan).clustered_peaks(2).peaks = sdtmp2.channel(chan).peaks;

    for clus = length(sdtmp.channel(chan).clustered_peaks):-1:1
        if isempty([sdtmp.channel(chan).clustered_peaks(clus).peaks.center_frequency])
            sdtmp.channel(chan).clustered_peaks(clus) = [];
            continue
        end
        if abs(mean([sdtmp.channel(chan).clustered_peaks(clus).peaks.center_frequency])-8)<1
            if pk1_seen
                i = i+1;
                sdtmp.channel(chan).clustered_peaks(clus) = [];
                continue
            end
            pk1_seen = 1;
            c = sdtmp.channel(chan).clustered_peaks(clus).peaks;
            if length([c.time]) ~= length(c)
                c(1:length(c)-length([c.time])) = [];
            end
            [~, iMat, iPk] = intersect(times(~isnan(pk1(1,:))),[c.time]);
            [~, iMat2, ~] = intersect(times,[c.time]);
            pk1err_cfr(chan,iMat) = [c(iPk).center_frequency];
            pk1err_amp(chan,iMat) = [c(iPk).amplitude];
            pk1err_std(chan,iMat) = [c(iPk).st_dev];
            pk1_det(chan,iMat2) = 1; 
            [~, iMat, iPk] = setxor(times(~isnan(pk1(1,:))),[c.time]);
            i = i+1;
            sdtmp.channel(chan).clustered_peaks(clus).peaks(iPk) = [];
        elseif abs(mean([sdtmp.channel(chan).clustered_peaks(clus).peaks.center_frequency])-16.5)<3
            if pk2_seen
                sdtmp.channel(chan).clustered_peaks(clus) = [];
                i = i+1;
                continue
            end
            pk2_seen = 1;
            c = sdtmp.channel(chan).clustered_peaks(clus).peaks;
            if length([c.time]) ~= length(c)
                c(1:length(c)-length([c.time])) = [];
            end
            [~, iMat, iPk] = intersect(times(~isnan(pk2(1,:))),[c.time]);
            [~, iMat2, ~] = intersect(times,[c.time]);
            pk2err_cfr(chan,iMat) = [c(iPk).center_frequency];
            pk2err_amp(chan,iMat) = [c(iPk).amplitude];
            pk2err_std(chan,iMat) = [c(iPk).st_dev];
            pk2_det(chan,iMat2) = 1; 
            [~, iMat, iPk] = setxor(times(~isnan(pk2(1,:))),[c.time]);
            i = i+1;
            sdtmp.channel(chan).clustered_peaks(clus).peaks(iPk) = [];
        end
    end
    pk1_seen = 0;
    pk2_seen = 0;
end
%%
% Formula for peak detection sensitivity and specificity
% For alpha (5.5-10.5 Hz)
alp_max_amp = pk1(2,:)./max(pk1(2,:))>=0.95; % maximum (>=95%) amplitude 
alp_max_amp_sens = sum(pk1_det(:,alp_max_amp),1);
mean_alp_max_amp_sens = mean(sum(pk1_det(:,alp_max_amp),1))./size(pk1_det,1); % alpha sensitivity
a = zeros(size(times)); a(pk1_exist) = 1;
alp_spec = sum(pk1_det(:,~a),1);
mean_alp_spec = 1-mean(alp_spec)./size(pk1_det,1); % alpha specificity
% For beta (13.5-20.5 Hz)
bet_max_amp = pk2(2,:)./max(pk2(2,:))>=0.95; % maximum (>=95%) amplitude 
bet_max_amp_sens = sum(pk2_det(:,bet_max_amp),1);
mean_bet_max_amp_sens = mean(sum(pk2_det(:,bet_max_amp),1))./size(pk2_det,1); % beta sensitivity
b = zeros(size(times)); b(pk2_exist) = 1;
bet_spec = sum(pk2_det(:,~a),1);
mean_bet_spec = 1-mean(bet_spec)./size(pk2_det,1); % beta specificity


%% Functions to help code run

function rel_heights = gen_peak_mat(len, fs, cf, amps, stds, tps, cfr, sdr, fooofFreqs, fooofTimes)

    freqs = linspace(0,fs./2,round(len.*fs./2));
    times = linspace(1./fs,len.*fs./fs,len.*fs);
    rel_heights = zeros(length(times),length(freqs));
    for peak = 1:length(cf)
        envMat = zeros(len*fs,1);
        % Generate centre frequency shift vector
        cfrMat = zeros(length(times),1); 
        cfrMat(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,2)+1:tps(peak,1)*fs+ ...
            tps(peak,2)*fs*cfr(peak,3))=linspace(0,cfr(peak,1),...
            round(cfr(peak,3)-cfr(peak,2),2)*tps(peak,2)*fs);
        cfrMat(tps(peak,1)*fs+tps(peak,2)*fs*cfr(peak,3):end) = cfr(peak,1);
        % Generate standard deviation shift vector
        sdrMat = zeros(size(times))';
        sdrMat(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,2)+1:tps(peak,1)*fs+ ...
            tps(peak,2)*fs*sdr(peak,3))=0;
        sdrMat(tps(peak,1)*fs+tps(peak,2)*fs*sdr(peak,3):end) = sdr(peak,1);
        envMat(tps(peak,1)*fs+1:(tps(peak,1)+tps(peak,2))*fs) = tukeywin(tps(peak,2).*fs,0.4);
        rh_tmp = ones(size(freqs)).*(envMat*amps(peak)).*...
            exp(-(ones(size(envMat)).*freqs-ones(size(freqs)).*...
            (cf(peak)+cfrMat)).^2./(2.*(ones(size(freqs)).*(stds(peak)+sdrMat).^2)));

        rel_heights = rel_heights +rh_tmp;
    end
    rh_tmp = zeros(length(freqs),length(fooofTimes));
    for f = 1:length(freqs)
        rh_tmp(f,:) = spline(times,rel_heights(:,f),fooofTimes);
    end
    rel_heights = zeros(length(fooofFreqs),length(fooofTimes));
    for t = 1:length(fooofTimes)
        rel_heights(:,t) = spline(freqs,rh_tmp(:,t),fooofFreqs);
    end
        
end

function ap_heights = gen_ap_mat1(slope, off, slr, fooofFreqs, fooofTimes)
    iTime1 = find(fooofTimes == slr(2));
    iTime2 = find(fooofTimes == slr(3));
    offMat = zeros(1,length(fooofTimes))+off(1); 
    offMat(iTime1:iTime2)=linspace(off(1), off(2), iTime2-iTime1+1);
    offMat(iTime2+1:end) = off(2);
    slrMat = zeros(1,length(fooofTimes))+slope; 
    slrMat(iTime1:iTime2)=linspace(slope, slr(1)+slope, iTime2-iTime1+1);
    slrMat(iTime2+1:end) = slope+slr(1);
    ap_heights = -log10(fooofFreqs)'*slrMat+ones(length(fooofFreqs),1)*offMat;
end

function draw_2575(times,measures,c)

x2 = [times fliplr(times)];
inBetween = [pchip(times,prctile(measures,25),times), fliplr(pchip(times,prctile(measures,75),times))];
f = fill(x2, inBetween, c ,'FaceAlpha',0.2);
f.EdgeColor = 'None';
f.FaceAlpha = 0.2;
hold on
plot(times, median(measures,'omitnan'), '-','Color',c)

end

function [h, u] = violin_plot(X, cl, offset,d,alp,rot,flipx)

[f, Xi, u] = ksdensity(X, linspace(min(X),max(X),1000));

Xii = [Xi(1) Xi Xi(end) fliplr(Xi)];
if flipx 
    fii = [0 f 0 zeros(1,length(f))].*-1./max(f).*d+offset;
else
    fii = [0 f 0 zeros(1,length(f))]./max(f).*d+offset;
end

if rot
    h{1} = fill(fii, Xii, cl); hold on
    h{1}.FaceAlpha = alp;
else
    h{1} = fill(Xii, fii, cl); hold on
    h{1}.FaceAlpha = alp;
end

set(h{1}, 'EdgeColor', [0.2 0.2 0.2]);

if rot 
    yl = get(gca, 'XLim');
    set(gca, 'XLim', [-yl(2) yl(2)]);
else
    yl = get(gca, 'YLim');
    set(gca, 'YLim', [-yl(2) yl(2)]);
end 
quartiles   = quantile(X, [0.25 0.75 0.5]);
iqr         = quartiles(2) - quartiles(1);
Xs          = sort(X);
whiskers(1) = min(Xs(Xs > (quartiles(1) - (1.5 * iqr))));
whiskers(2) = max(Xs(Xs < (quartiles(2) + (1.5 * iqr))));

plot([offset offset]-(-1+2.*flipx).*0.025, whiskers,'Color',cl,'LineWidth',0.2)
fill([offset offset offset-(-1+2.*flipx).*0.05 offset-(-1+2.*flipx).*0.05],quartiles([1 2 2 1]),[0.2 0.2 0.2],'EdgeColor','none')

scatter(offset-(-1+2.*flipx).*0.025,quartiles(3),8,[1 1 1],'Filled')

end

