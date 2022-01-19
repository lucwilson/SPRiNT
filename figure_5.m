%% Figure 5
% Focus on EC012 for now
work = whos('EC012_*'); % files are named EC012_1, EC012_2, etc.
work2 = whos('ec012ec*');
outstr = struct('velo',[]);

% Initialize PSDs
restTF12 = zeros(77,1);
ctr12 = 0;
movementTF12 = zeros(77,1);
ctm12 = 0;

% Initialize aperiodic parameter vectors
expsMat12 = [];
exprMat12 = [];
offsMat12 = [];
offrMat12 = [];

for t = 1:numel(work)
    sd_dat = eval(work(t).name);
    movdat = eval(work2(t).name);

    b = 1/39.06:1/39.06:length(movdat)/39.06;
    movdat(movdat(:,3)==-1,:)=nan;

    outstr(t).velo = sqrt((movdat(2:end,3)-movdat(1:end-1,3)).^2+(movdat(2:end,4)-movdat(1:end-1,4)).^2);
    outstr(t).pos = sqrt((movdat(:,3)-min(movdat(:,3))).^2+(movdat(:,4)-min(movdat(:,4))).^2);
    
    exp = zeros(length([sd_dat.tFOOOF.channel(1).data.time]),length(sd_dat.tFOOOF.channel));
    off = zeros(length([sd_dat.tFOOOF.channel(1).data.time]),length(sd_dat.tFOOOF.channel));
    error   = ones(length(sd_dat.Time),length(sd_dat.tFOOOF.channel));
    
    for d = 1:length(sd_dat.tFOOOF.channel)
        exp(:,d) = [sd_dat.tFOOOF.channel(d).aperiodics.exponent];
        off(:,d) = [sd_dat.tFOOOF.channel(d).aperiodics.offset];
        error(:,d) = sqrt([sd_dat.tFOOOF.channel(d).stats.MSE]);
    end
    exp = mean(exp,2);
    off = mean(off,2);
    
    c = b(~isnan(movdat(1:end-1,1)));
    outstr(t).exp = exp([sd_dat.tFOOOF.channel(1).data.time]>=ceil(c(1).*2)./2 &...
        [sd_dat.tFOOOF.channel(1).data.time]<=(floor(c(end).*2)./2)); 
    outstr(t).off = off([sd_dat.tFOOOF.channel(1).data.time]>=ceil(c(1).*2)./2 &...
        [sd_dat.tFOOOF.channel(1).data.time]<=(floor(c(end).*2)./2)); 
    outstr(t).error = error([sd_dat.tFOOOF.channel(1).data.time]>=ceil(c(1).*2)./2 &...
        [sd_dat.tFOOOF.channel(1).data.time]<=(floor(c(end).*2)./2));
    ts = [sd_dat.tFOOOF.channel(1).data.time];
    ts = ts(ts>=ceil(c(1).*2)./2 & ts<=(floor(c(end).*2)./2));
    tspos = ts(2:end)-(ts(2)-ts(1))./2;
    outstr(t).spdt = abs(diff(pchip(b,outstr(t).pos,tspos)));
    outstr(t).ts2 = ts(2:end-1);
    outstr(t).post = pchip(b,outstr(t).pos,outstr(t).ts2);
    
    outstr(t).velot = pchip(b(2:end),outstr(t).velo,ts);
    outstr(t).ts = ts;    
    
    spdthr = 6;

    indr = find(outstr(t).spdt>spdthr);
    inds = find(outstr(t).spdt<spdthr);

    % Add stationary points in middle of track to movement points
    indr = sort([inds(abs(outstr(t).post(inds)./max(outstr(t).post)-0.5) < 0.35) indr]);
    inds(abs(outstr(t).post(inds)./max(outstr(t).post)-0.5) < 0.35) = [];    

    % Add isolated movement points to stationary points
    for ii = length(indr):-1:1
        if length(indr) == ii
            if abs(indr(ii) - indr(ii-1)) > 1
                inds = [inds indr(ii)];
                indr(ii) = [];
            end
        elseif ii == 1
            if abs(indr(ii) - indr(ii+1)) > 1
                inds = [inds indr(ii)];
                indr(ii) = [];
            end
        else
              if (abs(indr(ii)-indr([ii-1 ii+1])) > 1)
                 inds = [inds indr(ii)];
                 indr(ii) = [];
             end
        end
    end
    % Locate rest (s) and movement (r) indices 
    inds = sort(inds);
    [~, ri, ~] = intersect(sd_dat.Time,outstr(t).ts2(indr));
    [~, si, ~] = intersect(sd_dat.Time,outstr(t).ts2(inds));
    % Combine all rest PSDs and movement PSDs
    ctr12 = ctr12 + length(si);
    restTF12 = restTF12 + squeeze(sum(mean(log10(sd_dat.TF(:,si,:)),1),2));
    ctm12 = ctm12 + length(ri);
    movementTF12 = movementTF12 + squeeze(sum(mean(log10(sd_dat.TF(:,ri,:)),1),2));    

    % Indexing requires shift one point forward
    inds = inds +1;
    indr = indr +1;
    % Set up boolean
    binds = zeros(length(outstr(t).ts2),1);
    binds(inds) = 1;
    bindr = zeros(length(outstr(t).ts2),1);
    bindr(indr) = 1;

    len = 8; % Minimum length of time bins before and after onset
    good_inds_rest = zeros(50,len.*2+1); % Rest onset indices
    ctr = 1;
    for ii = len:length(outstr(t).ts2) - len - 1
        if [bindr(ii-len+1:ii) ; binds((ii:ii+len)+1)]
            good_inds_rest(ctr,:) = ii-len+1:ii+len+1;
            ctr = ctr+1;
        end
    end
    good_inds_rest(ctr:end,:) = [];    
    
    good_inds_run = zeros(50,len.*2+1); % Movement onset indices
    ctr = 1;
    for ii = len:length(outstr(t).ts2) - len - 1
        if [binds(ii-len+1:ii) ; bindr((ii:ii+len)+1)]
            good_inds_run(ctr,:) = ii-len+1:ii+len+1;
            ctr = ctr+1;
        end
    end
    good_inds_run(ctr:end,:) = [];
    % Save indices, matrices of aperiodic parameters
    outstr(t).good_inds_rest = good_inds_rest;
    outstr(t).good_inds_run = good_inds_run;
    outstr(t).indr = indr;
    outstr(t).inds = inds;
    expsMat12 = [expsMat12; outstr(t).exp(inds)];
    exprMat12 = [exprMat12; outstr(t).exp(indr)];
    offsMat12 = [offsMat12; outstr(t).off(inds)];
    offrMat12 = [offrMat12; outstr(t).off(indr)];
end
% Divide by number of spectra in each run and movement block.
restTF12 = restTF12./ctr12;
movementTF12 = movementTF12./ctm12;
% Extract aperiodic parameters throughout epochs
good12_rest = zeros(400,len.*2+1);
good12_run = zeros(400,len.*2+1);
ctr = 1; % For rest
for h = 1:9
    for q = 1:size(outstr(h).good_inds_rest,1)
        good12_rest(ctr,:) = outstr(h).exp(outstr(h).good_inds_rest(q,:));
        ctr = ctr+1;
    end
end
good12_rest(ctr:end,:) = [];
ctr = 1; % For movement
for h = 1:9
    for q = 1:size(outstr(h).good_inds_run,1)
        good12_run(ctr,:) = outstr(h).exp(outstr(h).good_inds_run(q,:));
        ctr = ctr+1;
    end
end
good12_run(ctr:end,:) = [];
outstr12 = outstr;



% Now focus on EC013
work = whos('EC013_*'); % files are named EC013_1, EC013_2, etc.
work2 = whos('ec013ec*');
outstr = struct('velo',[]);

% Initialize PSDs
restTF13 = zeros(77,1);
ctr13 = 0;
movementTF13 = zeros(77,1);
ctm13 = 0;

% Initialize aperiodic parameter vectors
expsMat13 = [];
exprMat13 = [];
offsMat13 = [];
offrMat13 = [];

for t = 1:numel(work)
    sd_dat = eval(work(t).name);
    movdat = eval(work2(t).name);

    b = 1/39.06:1/39.06:length(movdat)/39.06;
    movdat(movdat(:,3)==-1,:)=nan;
    a = fitlm(movdat(:,3),movdat(:,4));

    outstr(t).velo = sqrt((movdat(2:end,3)-movdat(1:end-1,3)).^2+(movdat(2:end,4)-movdat(1:end-1,4)).^2);
    
    if a.Coefficients{2,1} > 0
        outstr(t).pos = sqrt((movdat(:,3)-min(movdat(:,3))).^2+(movdat(:,4)-min(movdat(:,4))).^2);
    else
        outstr(t).pos = sqrt((movdat(:,3)-max(movdat(:,3))).^2+(movdat(:,4)-min(movdat(:,4))).^2);
    end
    exp = zeros(length([sd_dat.tFOOOF.channel(1).data.time]),length(sd_dat.tFOOOF.channel));
    off = zeros(length([sd_dat.tFOOOF.channel(1).data.time]),length(sd_dat.tFOOOF.channel));
    error   = ones(length(sd_dat.Time),length(sd_dat.tFOOOF.channel));
    
    for d = 1:length(sd_dat.tFOOOF.channel)
        exp(:,d) = [sd_dat.tFOOOF.channel(d).aperiodics.exponent];
        off(:,d) = [sd_dat.tFOOOF.channel(d).aperiodics.offset];
        error(:,d) = sqrt([sd_dat.tFOOOF.channel(d).stats.MSE]);
    end
    exp = mean(exp,2);
    off = mean(off,2);
    
    c = b(~isnan(movdat(1:end-1,1)));
    outstr(t).exp = exp([sd_dat.tFOOOF.channel(1).data.time]>=ceil(c(1).*2)./2 &...
        [sd_dat.tFOOOF.channel(1).data.time]<=(floor(c(end).*2)./2)); 
    outstr(t).off = off([sd_dat.tFOOOF.channel(1).data.time]>=ceil(c(1).*2)./2 &...
        [sd_dat.tFOOOF.channel(1).data.time]<=(floor(c(end).*2)./2)); 
    outstr(t).error = error([sd_dat.tFOOOF.channel(1).data.time]>=ceil(c(1).*2)./2 &...
        [sd_dat.tFOOOF.channel(1).data.time]<=(floor(c(end).*2)./2));
    ts = [sd_dat.tFOOOF.channel(1).data.time];
    ts = ts(ts>=ceil(c(1).*2)./2 & ts<=(floor(c(end).*2)./2));
    tspos = ts(2:end)-(ts(2)-ts(1))./2;
    outstr(t).spdt = abs(diff(pchip(b,outstr(t).pos,tspos)));
    outstr(t).ts2 = ts(2:end-1);
    outstr(t).post = pchip(b,outstr(t).pos,outstr(t).ts2);
    
    outstr(t).velot = pchip(b(2:end),outstr(t).velo,ts);
    outstr(t).ts = ts;    
    
    spdthr = 6;

    indr = find(outstr(t).spdt>spdthr);
    inds = find(outstr(t).spdt<spdthr);

    % Add stationary points in middle of track to movement points
    indr = sort([inds(abs(outstr(t).post(inds)./max(outstr(t).post)-0.5) < 0.35) indr]);
    inds(abs(outstr(t).post(inds)./max(outstr(t).post)-0.5) < 0.35) = [];    

    % Add isolated movement points to stationary points
    for ii = length(indr):-1:1
        if length(indr) == ii
            if abs(indr(ii) - indr(ii-1)) > 1
                inds = [inds indr(ii)];
                indr(ii) = [];
            end
        elseif ii == 1
            if abs(indr(ii) - indr(ii+1)) > 1
                inds = [inds indr(ii)];
                indr(ii) = [];
            end
        else
              if (abs(indr(ii)-indr([ii-1 ii+1])) > 1)
                 inds = [inds indr(ii)];
                 indr(ii) = [];
             end
        end
    end
    % Locate rest (s) and movement (r) indices 
    inds = sort(inds);
    [~, ri, ~] = intersect(sd_dat.Time,outstr(t).ts2(indr));
    [~, si, ~] = intersect(sd_dat.Time,outstr(t).ts2(inds));
    % Combine all rest PSDs and movement PSDs
    ctr13 = ctr13 + length(si);
    restTF13 = restTF13 + squeeze(sum(mean(log10(sd_dat.TF(:,si,:)),1),2));
    ctm13 = ctm13 + length(ri);
    movementTF13 = movementTF13 + squeeze(sum(mean(log10(sd_dat.TF(:,ri,:)),1),2));    

    % Indexing requires shift one point forward
    inds = inds +1;
    indr = indr +1;
    % Set up boolean
    binds = zeros(length(outstr(t).ts2),1);
    binds(inds) = 1;
    bindr = zeros(length(outstr(t).ts2),1);
    bindr(indr) = 1;

    len = 8; % Minimum length of time bins before and after onset
    good_inds_rest = zeros(50,len.*2+1); % Rest onset indices
    ctr = 1;
    for ii = len:length(outstr(t).ts2) - len - 1
        if [bindr(ii-len+1:ii) ; binds((ii:ii+len)+1)]
            good_inds_rest(ctr,:) = ii-len+1:ii+len+1;
            ctr = ctr+1;
        end
    end
    good_inds_rest(ctr:end,:) = [];    
    
    good_inds_run = zeros(50,len.*2+1); % Movement onset indices
    ctr = 1;
    for ii = len:length(outstr(t).ts2) - len - 1
        if [binds(ii-len+1:ii) ; bindr((ii:ii+len)+1)]
            good_inds_run(ctr,:) = ii-len+1:ii+len+1;
            ctr = ctr+1;
        end
    end
    good_inds_run(ctr:end,:) = [];
    % Save indices, matrices of aperiodic parameters
    outstr(t).good_inds_rest = good_inds_rest;
    outstr(t).good_inds_run = good_inds_run;
    outstr(t).indr = indr;
    outstr(t).inds = inds;
    expsMat13 = [expsMat13; outstr(t).exp(inds)];
    exprMat13 = [exprMat13; outstr(t).exp(indr)];
    offsMat13 = [offsMat13; outstr(t).off(inds)];
    offrMat13 = [offrMat13; outstr(t).off(indr)];
end
% Divide by number of spectra in each run and movement block.
restTF13 = restTF13./ctr13;
movementTF13 = movementTF13./ctm13;
% Extract aperiodic parameters throughout epochs
good13_rest = zeros(400,len.*2+1);
good13_run = zeros(400,len.*2+1);
ctr = 1; % For rest
for h = 1:9
    for q = 1:size(outstr(h).good_inds_rest,1)
        good13_rest(ctr,:) = outstr(h).exp(outstr(h).good_inds_rest(q,:));
        ctr = ctr+1;
    end
end
good13_rest(ctr:end,:) = [];
ctr = 1; % For movement
for h = 1:9
    for q = 1:size(outstr(h).good_inds_run,1)
        good13_run(ctr,:) = outstr(h).exp(outstr(h).good_inds_run(q,:));
        ctr = ctr+1;
    end
end
good13_run(ctr:end,:) = [];
outstr13 = outstr;


% Figure 5
figure('Position',[200 200 900 600])

c1 = [1 50 150]./255; 
c2 = [239 161 7]./255;

% Rat 1
subplot(4,3,4), hold on
plot(restTF12,2:0.5:40,'Color',c2,'LineWidth',1.2)
plot(movementTF12,2:0.5:40,'Color',c1,'LineWidth',1.2)
ylabel('Frequency (Hz)')
xlabel('Power (a.u.)')
set(gca,'YAxisLocation','right','xdir','reverse');
ylim([2 40])
yticks(10:10:30)
xlim([-4 -0.7])

subplot(4,3,5), hold on
imagesc(EC012_6.Time(217:337),2:0.5:40,squeeze(mean(log10(EC012_6.tFOOOF.tFOOOF_models(:,217:337,:)),1))')
set(gca,'ydir','normal');
ylim([2 40])
xlim(EC012_6.Time([217 337]))
colorbar;
caxis([-4 -0.7])
ylabel('Frequency (Hz)')
xlabel('Time (s)')
xticks(100:20:300)
yticks(10:10:30)

subplot(4,3,2), hold on
spdthr = 5;
bi = zeros(length(outstr12(6).ts2),1);
bi(outstr12(6).indr) = 6;
draw_shade(outstr12(6).ts2,bi,spdthr,1,c1,0.3)
draw_shade(outstr12(6).ts2,bi,spdthr,0,c2,0.3)
plot(1/39.06:1/39.06:length(outstr12(6).pos)/39.06,outstr12(6).pos./max(outstr12(6).pos),'k');
xlim(EC012_6.Time([217 337]))
xticks(100:20:300)
xticklabels([])
ylabel('Position (cm)')
xlabel('Time (s)')
yticks(0:0.5:1)
yticklabels(0:125:250)

subplot(4,3,6), hold on
plot(-4:0.5:4,mean(good12_rest),'Color',c2)
plot(-4:0.5:4,mean(good12_run),'Color',c1)
draw_95ci(-4:0.5:4,good12_rest)
draw_95ci(-4:0.5:4,good12_run)
plot(-4:0.5:4,mean(good12_rest),'Color',c2,'LineWidth',1.2)
plot(-4:0.5:4,mean(good12_run),'Color',c1,'LineWidth',1.2)
legend('0 = rest onset','0 = movement onset')
xlabel('Time (s)')
ylabel('Exponent (Hz^{-1})')
ylim([0.95 2.05])
yticks(1:0.25:2)


% Rat 2

subplot(4,3,10), hold on
plot(restTF13,2:0.5:40,'Color',c2,'LineWidth',1.2)
plot(movementTF13,2:0.5:40,'Color',c1,'LineWidth',1.2)
ylabel('Frequency (Hz)')
xlabel('Power (a.u.)')
set(gca,'YAxisLocation','right','xdir','reverse');
yticks(10:10:30)
xlim([-4 -0.7])
ylim([2 40])

subplot(4,3,11), hold on
imagesc(EC013_7.Time(217:337),2:0.5:40,squeeze(mean(log10(EC013_7.tFOOOF.tFOOOF_models(:,217:337,:)),1))')
set(gca,'ydir','normal');
caxis([-4 -0.7])
colorbar;
xlim(EC013_7.Time([217 337]))
ylim([2 40])
yticks(10:10:30)
ylabel('Frequency (Hz)')
xlabel('Time (s)')
xticks(100:20:300)

subplot(4,3,8), hold on
spdthr = 5;
bi = zeros(length(outstr13(7).ts2),1);
bi(outstr13(7).indr) = 6;
draw_shade(outstr13(7).ts2,bi,spdthr,1,c1,0.3)
draw_shade(outstr13(7).ts2,bi,spdthr,0,c2,0.3)
plot(1/39.06:1/39.06:length(outstr13(7).pos)/39.06,outstr13(7).pos./max(outstr13(7).pos),'k');
xlim(EC013_7.Time([217 337]))
xticks(100:20:300)
% xticklabels([])
ylabel('Position (cm)')
xlabel('Time (s)')
yticks(0:0.5:1)
yticklabels(0:125:250)

subplot(4,3,12), hold on
plot(-4:0.5:4,mean(good13_rest),'Color',c2)
plot(-4:0.5:4,mean(good13_run),'Color',c1)
draw_95ci(-4:0.5:4,good13_rest)
draw_95ci(-4:0.5:4,good13_run)
plot(-4:0.5:4,mean(good13_rest),'Color',c2,'LineWidth',1.2)
plot(-4:0.5:4,mean(good13_run),'Color',c1,'LineWidth',1.2)
legend('0 = rest onset','0 = movement onset')
xlabel('Time (s)')
ylabel('Exponent (Hz^{-1})')
ylim([0.95 2.05])
yticks(1:0.25:2)


%% Support functions
function draw_shade(times, speed, thresh, above, cl, alp)
    sRate = times(2)-times(1);
    if above
        indgood = speed > thresh;
    else
        indgood = speed <= thresh;
    end
    drawing = 0;
    ts = [0 0];
    vals = [0.9 0.9 1 1];
    for i = 1:length(times)
        if indgood(i)
            if ~drawing
                ts(1) = times(i)-sRate./2;
                drawing = 1;
            else
                continue
            end
        else
            drawing = 0;
            if i > 1
                if indgood(i-1)
                    ts(2) = times(i-1)+sRate./2;
                    drt = [ts fliplr(ts)];
                    a = fill(drt,vals,cl);
                    a.FaceAlpha = alp;
                    a.EdgeColor = 'None';
                end
            end
        end
    end
end

function draw_95ci(times,measures)
x2 = [times fliplr(times)];
tval = tinv(0.025,size(measures,1));
inBetween = [mean(measures)-tval.*std(measures)./sqrt(size(measures,1)),...
    fliplr(mean(measures)+tval.*std(measures)./sqrt(size(measures,1)))];
f = fill(x2, inBetween, [205 205 205]./255,'FaceAlpha',0.4);
f.EdgeColor = 'None';
f.FaceAlpha = 0.5;
plot(times, mean(measures,'omitnan'), '-b')
plot(times, mean(measures)-tval.*std(measures)./sqrt(size(measures,1)), '-k')
plot(times, mean(measures)+tval.*std(measures)./sqrt(size(measures,1)), '-k')
end