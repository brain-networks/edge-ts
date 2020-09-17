clear all
close all
clc

% add helper functions to path
addpath(genpath('fcn'));

% load regional time series
load data/ts.mat

% z-score time series
ts = zscore(ts);

% get dimensions
[ntime,nnodes] = size(ts);

% calculate number of edges
nedges = nnodes*(nnodes - 1)/2;

% indices of unique edges (upper triangle)
[u,v] = find(triu(ones(nnodes),1));
idx = (v - 1)*nnodes + u;

%% calculate static fc
fc = corr(ts);

%% generate edge time series
ets = ts(:,u).*ts(:,v);

%% plot average time averaged edge time series against upper triangle fc
fc_upper_triangle = fc(idx);
figure,plot(fc_upper_triangle,mean(ets),'ko')
xlabel('upper triangle fc matrix');
ylabel('time average of edge time series');
xlim([-1,1]); ylim([-1,1]);
text(0.8,-0.8,sprintf('r = %.2f',corr(fc_upper_triangle,mean(ets)')));

%% calculate event amplitude

% calculate co-fluctuation amplitude at each frame
rms = sum(ets.^2,2).^0.5;

% fraction of high- and low-amplitude frames to retain
frackeep = 0.1;
nkeep = round(ntime*frackeep);

% sort co-fluctuation amplitude
[~,idxsort] = sort(rms,'descend');

% estimate fc using just high-amplitude frames
fctop = corr(ts(idxsort(1:nkeep),:));

% do the same using just low-amplitude
fcbot = corr(ts(idxsort(end - nkeep + 1:end),:));


figure('position',[560,528,700,420])
subplot(3,2,2); plot(rms);
xlabel('frames'); ylabel('rss');

subplot(3,2,1); imagesc(fc,[-1,1]); axis square;
title('fc all time points');

subplot(3,2,3); imagesc(fctop,[-1,1]); axis square;
title('fc high-amplitude time points');

subplot(3,2,4); scatter(fctop(idx),fc(idx),'k.'); xlim([-1,1]); ylim([-1,1]);
text(-1,1,sprintf('r = %.2f',corr(fctop(idx),fc(idx))));

subplot(3,2,5); imagesc(fcbot,[-1,1]); axis square;
title('fc low-amplitude time points');

subplot(3,2,6); scatter(fcbot(idx),fc(idx),'k.'); xlim([-1,1]); ylim([-1,1]);
text(-1,1,sprintf('r = %.2f',corr(fcbot(idx),fc(idx))));

%% calculate similarity of time-averaged fc with high-/low-amplitude frames
rho(1) = corr(fctop(idx),fc(idx));
rho(2) = corr(fcbot(idx),fc(idx));

%% calculate modularity of high-/low-amplitude fc
numiter = 100;
qtop = zeros(numiter,1); qbot = qtop;
for iter = 1:numiter
    [~,qtop(iter)] = community_louvain(fctop,[],[],'negative_asym');
    [~,qbot(iter)] = community_louvain(fcbot,[],[],'negative_asym');
end
figure;
f = fcn_boxpts([qtop,qbot],[],jet(2));
set(gca,'xtick',1:2,'xticklabel',{'high-amp','low-amp'});
ylabel('modularity, q');