function evaluate_decoder(learning_dataset, output_dataset)
lrn = load_runtime_data(learning_dataset);
opt = load_runtime_data(output_dataset);
figure; set(gcf, 'color', 'w');
% plot the two curves and numerical RMSE
subplot(2,1,1);
plot(lrn.sim.indata.data(:,1),lrn.sim.indata.data(:,2),'.r'); hold on;
plot(opt.sim.indata.data(:,1),opt.sim.indata.data(:,2),'.g'); hold on; box off;
rmse_opt = sqrt(sum((lrn.sim.indata.data(:,2) - opt.sim.indata.data(:,2)).^2)/numel(lrn.sim.indata.data(:,2)));
legend('Learned function',sprintf('Optimizer based decoder recovery - RMSE : %f', rmse_opt));
% plot the input data and the wrapping manifold depicting decoded error
subplot(2,1,2);
errorbar(lrn.sim.indata.data(:,1), lrn.sim.indata.data(:,2), lrn.sim.indata.data(:,2) - opt.sim.indata.data(:,2), '.g')
box off; hold on; plot(lrn.sim.indata.data(:,1),lrn.sim.indata.data(:,2),'.r');
end