function [f,hdl,pts] = fcn_boxpts(x,lab,cmap,samenum,w)
% clear all
% close all
% clc

% load boxdata.mat
% x = mu;

if ~exist('samenum','var') | isempty(samenum)
    samenum = 0;
end
h = hist(lab,1:max(lab));

if ~exist('w','var') | isempty(samenum)
    w = 1.5;
end

sz = size(x);
if sz(2) > 1
    lab = repmat(1:sz(2),sz(1),1);
    x = x(:);
    lab = lab(:);
end

if isempty(cmap)
    cmap = jet(max(lab));
end

iqr = [25,50,75];
wid = 0.3;
cmapdark = cmap - 0.25;
cmapdark(cmapdark < 0) = 0;

u = unique(lab);

f = figure(gcf);
ax = axes;
hold(ax,'on');

for j = 1:length(u)
    idx = lab == u(j);
    vals = x(idx);
    prct = prctile(vals,iqr);
    
    dr = prct(3) - prct(1);
    upper = prct(3) + w*dr;
    lower = prct(1) - w*dr;
    
    whisk1 = min(vals(vals >= lower));
    whisk2 = max(vals(vals <= upper));
    
    if samenum ~= 0
        r = randperm(length(vals),samenum);
        vals = vals(r);
    end
    
    xvals = j + linspace(-wid*0.5,wid*0.5,length(vals));
    s = randperm(length(xvals));
    xvals = xvals(s);
    
    hdl.pointshandle{j} = scatter(xvals,vals,5,'filled');
    set(hdl.pointshandle{j},'markerfacecolor',cmap(j,:),'markerfacealpha',0.5);
    
    pts{j}=[xvals(:),vals(:)];
    
    xvals = j + [-wid,+wid,+wid,-wid,-wid];
    yvals = [prct(1),prct(1),prct(3),prct(3),prct(1)];
    hdl.boxhandle{j} = plot(xvals,yvals,'color','k','linewidth',0.75);
    hdl.medianhandle{j} = plot(j + [-wid,wid],prct(2)*ones(1,2),'color','k','linewidth',0.75);
    
    
    xvals = [j*ones(1,2), nan, j*ones(1,2)];
    yvals = [prct(1),whisk1,nan,prct(3),whisk2];
    hdl.whiskerhandle{j} = plot(xvals,yvals,'color','k','linewidth',0.75);
    
end
set(gca,...
    'xlim',[-2*wid + 1,length(u) + 2*wid],...
    'ylim',[min(x) - 0.05*range(x),max(x) + 0.05*range(x)])

