    % load the data from the file
function visualize_runtime_nst_thesis(filein)
rdata = load_runtime_data(filein);
close all;
% plot runtime and learning parameters of the network
figure;
set(gcf, 'color', 'white');
subplot(2,2,1);
plot(rdata.sim.alpha, '.k'); xlabel('in learning epochs'); ylabel('alpha');
grid off; box off; title(sprintf('Input learning rate, adapt %d epochs', rdata.sim.tf_lrn_in));
subplot(2,2,2);
plot(rdata.sim.sigma, '.k'); xlabel('in learning epochs'); ylabel('sigma');
grid off; box off; title(sprintf('Neighborhood size, adapt %d epochs', rdata.sim.tf_lrn_in));
subplot(2,2,3);
plot(rdata.sim.eta, '.k'); xlabel('cross learning epochs'); ylabel('eta');
grid off; box off; title(sprintf('Activity decay factor, adapt %d epochs', rdata.sim.tf_lrn_cross));
subplot(2,2,4);
plot(rdata.sim.xi, '.k'); xlabel('cross learning epochs'); ylabel('xi');
grid off; box off; title(sprintf('Cross-modal learning rate, adapt %d epochs', rdata.sim.tf_lrn_cross));

% figure for nst thesis
idcolor = 'rbgmyk';
for pidx = 2:rdata.sim.net.nsize
    % input data and analysis
    figure; set(gcf, 'color','w');
    subplot(8, 8, [3 14]);
    plot(rdata.sim.indata.data(:, 1), 'r'); box off; set(gca,'YTick',[], 'XTick',[],'linewidth',3);
    subplot(8, 8, [17,42]);
    [c, b]=hist(rdata.sim.indata.data(:, pidx), 50);set(gca,'YTick',[], 'XTick',[], 'xdir', 'n','linewidth',3);
    h = barh(b, c, 'hist');box off;
    set(get(h, 'Parent'), 'xdir', 'r');  set(gca,'visible','off');
    subplot(8, 8, [19 46]);
    plot(rdata.sim.indata.data(:, 1), rdata.sim.indata.data(:, pidx), '.k'); box off;set(gca,'YTick',[], 'XTick',[],'linewidth',3);
    subplot(8, 8, [23 48]);
    plot(rdata.sim.indata.data(:, pidx), sprintf('%s',idcolor(pidx))); box off; set(gca,'YTick',[], 'XTick',[],'linewidth',3);
    view([90 90]);
    subplot(8, 8, [51 62]);
    hist(rdata.sim.indata.data(:, 1), 50);set(gca,'visible','off');
end

idcolor = 'bgmyk';
for pidx = 2:rdata.sim.net.nsize
    % learned data and analysis
    figure; set(gcf, 'color','w');
    subplot(8, 8, [3 14]);
    plot(rdata.sim.net.pops(1).s, '.r'); box off; set(gca,'YTick',[], 'XTick',[],'linewidth',3);
    subplot(8, 8, [17,42]);
    [c, b]=hist(rdata.sim.net.pops(pidx).Winput, 50);
    h = barh(b, c, 'hist'); box off;
    set(get(h, 'Parent'), 'xdir', 'r');  set(gca,'visible','off');
    subplot(8, 8, [19 46]);
    imagesc((rdata.sim.net.pops(pidx).Wcross), [0, 1]); box off; set(gca,'YTick',[], 'XTick',[]);
    subplot(8, 8, [23 48]);
    plot(rdata.sim.net.pops(pidx).s, sprintf('.%s',idcolor(pidx))); box off; view([90 90]);set(gca,'YTick',[], 'XTick',[],'linewidth',3);
    subplot(8, 8, [51 62]);
    hist(rdata.sim.net.pops(1).Winput, 50); set(gca,'visible','off');
end
end