% load the data from the file
function visualize_runtime_nst_thesis(filein, dec_filein)
rdata = load_runtime_data(filein);
trdata = load_runtime_data(dec_filein);
% %close all;
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
    plot(rdata.sim.indata.data(:, 1), rdata.sim.indata.data(:, pidx), '.k'); box off;%set(gca,'YTick',[], 'XTick',[],'linewidth',3);
    subplot(8, 8, [23 48]);
    plot(rdata.sim.indata.data(:, pidx), sprintf('%s',idcolor(pidx))); box off; set(gca,'YTick',[], 'XTick',[],'linewidth',3);
    view([90 90]);
    subplot(8, 8, [51 62]);
    hist(rdata.sim.indata.data(:, 1), 50);set(gca,'visible','off');
end

idcolor = 'rbgmyk';
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
    imagesc(rot90(rdata.sim.net.pops(pidx).Wcross), [0, 1]); box off; set(gca,'YTick',[], 'XTick',[]);
    subplot(8, 8, [23 48]);
    plot(rdata.sim.net.pops(pidx).s, sprintf('.%s',idcolor(pidx))); box off; view([90 90]);set(gca,'YTick',[], 'XTick',[],'linewidth',3);
    subplot(8, 8, [51 62]);
    hist(rdata.sim.net.pops(1).Winput, 50); set(gca,'visible','off');
end
% individual plots for analysis: input data and decoded data
figure; set(gcf, 'color','w');
for pidx=2:rdata.sim.net.nsize
    subplot(1, rdata.sim.net.nsize-1, pidx-1);
    plot(rdata.sim.indata.data(:, 1), rdata.sim.indata.data(:, pidx), '.k', 'linewidth',3); box off;
    hold on;
    plot(trdata.sim.indata.data(:, 1), trdata.sim.indata.data(:, pidx), '.c', 'linewidth',3); box off;
end
% individual plots for analysis: learned coractivation pattern
figure; set(gcf, 'color','w');
for pidx=2:rdata.sim.net.nsize
    subplot(1, rdata.sim.net.nsize-1, pidx-1);
    switch(pidx)
        case 2
            imagesc(rot90(rdata.sim.net.pops(pidx).Wcross), [0, 1]);
        case 3
            imagesc(rot90(rdata.sim.net.pops(pidx).Wcross'), [0, 1]);
        case 4
            imagesc((rdata.sim.net.pops(pidx).Wcross), [0, 1]);
    end
    %box off; set(gca,'YTick',[], 'XTick',[]);
end
% separate plots for thesis analysis on multimodal extensions of the model
figure; set(gcf, 'color','w');
plot(rdata.sim.indata.data(:, 1), rdata.sim.indata.data(:, 2), '.k', 'linewidth',3); box off; 
hold on;
plot(trdata.sim.indata.data(:, 1), trdata.sim.indata.data(:, 2), '.c', 'linewidth',3); box off;
set(gca,'YTick',[], 'XTick',[]);
figure; set(gcf, 'color','w');
imagesc((rdata.sim.net.pops(1).Wcross'), [0, 1]); box off; set(gca,'YTick',[], 'XTick',[]);
end