    % load the data from the file
function visualize_runtime(filein)
rdata = load_runtime_data(filein);
close all;
% plot runtime and learning parameters of the network
figure(1);
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

% plot runtime and learning parameters of the network
idcolor = 'rbgmyk';
for pidx = 1:rdata.sim.net.nsize
    figure;
    set(gcf, 'color', 'white');
    subplot(3, 1, 1);
    plot(rdata.sim.net.pops(pidx).Winput, sprintf('.%s',idcolor(pidx))); xlabel('cross learning epochs');
    ylabel(sprintf('Winput - pop %d', pidx));
    grid off; box off; title(sprintf('Input weight vector, adapt %d epochs', rdata.sim.tf_lrn_in));
    subplot(3, 1, 2);
    imagesc(rot90(rdata.sim.net.pops(pidx).Wcross), [0, 1]); box off; colorbar;
    xlabel('cross learning epochs'); ylabel(sprintf('Wcross - pop %d', pidx));
    grid off; box off; title(sprintf('Cross weight vector, adapt %d epochs', rdata.sim.tf_lrn_cross));
    subplot(3, 1, 3);
    plot(rdata.sim.net.pops(pidx).s, sprintf('.%s',idcolor(pidx)));
    xlabel('cross learning epochs'); ylabel(sprintf('tuning curve size, s - pop %d', pidx));
    grid off; box off; title(sprintf('Tuning curves values, adapt %d epochs', rdata.sim.tf_lrn_cross));
end

% plot input data
RANGE = 1;
for vpidx =2:rdata.sim.net.nsize
    figure; set(gcf, 'color', 'white');
    plot(rdata.sim.indata.data(:, 1), rdata.sim.indata.data(:, vpidx), '.k'); xlabel('samples');
    grid off; box off; title(sprintf('Input relation population 1 - population %d', vpidx));
end
% individual map analysis
for ppidx = 1:rdata.sim.indata.npop
    figure;
    set(gcf, 'color', 'white');
    subplot(4, 1, 1);
    plot(rdata.sim.indata.data(:, ppidx), sprintf('%s',idcolor(ppidx))); xlabel('samples');
    grid off; box off; title(sprintf('Input data - pop %d', ppidx));
    subplot(4, 1, 2);
    histogram(rdata.sim.indata.data(:, ppidx), 50, 'FaceColor',sprintf('%s',idcolor(ppidx)));
    grid off; box off; title(sprintf('Input data distribution - pop %d', ppidx));
    hndl = subplot(4, 1, 3);
    % compute the tuning curve of the current neuron in the population
    % the equally spaced mean values
    x = linspace(-RANGE, RANGE, rdata.sim.indata.popsize);
    % for each neuron in the current population compute the receptive field
    for idx = 1:rdata.sim.indata.popsize
        % extract the preferred values (wight vector) of each neuron
        v_pref = rdata.sim.net.pops(ppidx).Winput(idx);
        fx = exp(-(x - v_pref).^2/(2*rdata.sim.net.pops(ppidx).s(idx)^2));
        plot(1:rdata.sim.indata.popsize, fx, 'LineWidth', 3); hold all;
    end
    rdata.sim.net.pops(ppidx).Winput = sort(rdata.sim.net.pops(ppidx).Winput); box off;
    ax1_pos = get(hndl, 'Position'); set(hndl, 'XTick', []); set(hndl, 'XColor','w');
    ax2 = axes('Position',ax1_pos,'XAxisLocation','bottom','Color','none','LineWidth', 3);
    set(hndl, 'YTick', []); set(hndl, 'YColor','w');
    set(ax2, 'XTick', rdata.sim.net.pops(ppidx).Winput); set(ax2, 'XTickLabel', []);
    set(ax2, 'XLim', [ min(rdata.sim.net.pops(ppidx).Winput), max(rdata.sim.net.pops(ppidx).Winput)]);
    xlabel('neuron preferred values'); title('Learned tuning curves shapes');
    % the density of the tuning curves (density function) - should increase
    % with the increase of the distribution of sensory data (directly proportional with the prior, p(s))
    % stimuli associated with the peaks of the tuning curves
    subplot(4, 1, 4);
    histogram(rdata.sim.net.pops(ppidx).Winput, 50, 'FaceColor',sprintf('%s',idcolor(ppidx))); box off;title('# of allocated neurons for a value');
    xlabel('input value range');
end

% DISPLAY TUNING CURVES FOR ONLY SOME NEURONS
% individual map analysis
for ppidx = 1:rdata.sim.indata.npop
figure;
hndl = subplot(1,1,1);
v_pref = sort(rdata.sim.net.pops(ppidx).Winput);
% for each neuron in the current population compute the receptive field
% select some tuning curves to plot
pref = [1, 6, 13, 40, 45, 85, 90, 99];
for idx = 1:length(pref)
    idx_pref = pref(idx);
    % extract the preferred values (weight vector) of each neuron
    fx = exp(-(x - v_pref(idx_pref)).^2/(2*rdata.sim.net.pops(ppidx).s(idx_pref)^2));
    plot(1:rdata.sim.indata.popsize, fx,'LineWidth', 3); hold all;
end
ax1_pos = get(hndl, 'Position'); set(hndl, 'XTick', []); set(hndl, 'XColor','w');
ax2 = axes('Position',ax1_pos,'XAxisLocation','bottom','Color','none','LineWidth', 3);
set(hndl, 'YTick', []); set(hndl, 'YColor','w');
v_pref_idx = zeros(1, length(pref));
for idx = 1:length(pref)
    v_pref_idx(idx) = v_pref(pref(idx)+1);
end
set(ax2, 'XTick', v_pref_idx); set(ax2, 'XTickLabel', v_pref_idx);
set(ax2, 'XLim', [min(x), max(x)]);
set(ax2, 'XTickLabelRotation', 90);
xlabel(ax2, 'preferred value');
ylabel('learned tuning curves shapes');
ax3 = axes('Position',ax2.Position,...
    'XAxisLocation','top',...
    'Color','none','LineWidth', 0.01);
set(ax3, 'XTick', v_pref_idx); 
set(ax3, 'Ticklength', [0 0]); 
set(ax3, 'XTickLabel', pref);
set(ax3, 'XLim', [min(x), max(x)]);
xlabel(ax3, 'neuron index');
end
end