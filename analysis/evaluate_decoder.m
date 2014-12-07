close all;
lrn = load_runtime_data('../learning_data');
opt = load_runtime_data('../output_decoded_optimizer');
naive = load_runtime_data('../output_decoded_naive');
figure; set(gcf, 'color', 'w');
plot(lrn.sim.indata.data(:,1),lrn.sim.indata.data(:,2),'.r'); hold on;
plot(opt.sim.indata.data(:,1),opt.sim.indata.data(:,2),'.g'); hold on;
plot(naive.sim.indata.data(:,1),naive.sim.indata.data(:,2),'.b'); hold off; box off;
legend('Learned function','Optimizer based decoder recovery','Naive decoder recovery');
figure; set(gcf, 'color', 'w');
rmse_opt = sqrt(sum((lrn.sim.indata.data(:,2) - opt.sim.indata.data(:,2)).^2)/numel(lrn.sim.indata.data(:,2)));
rmse_naive = sqrt(sum((lrn.sim.indata.data(:,2) - naive.sim.indata.data(:,2)).^2)/numel(lrn.sim.indata.data(:,2)));
rmse_opt_v = zeros(1, length(lrn.sim.indata.data(:,2)));
rmse_naive_v = zeros(1, length(lrn.sim.indata.data(:,2)));
for i=1:length(lrn.sim.indata.data(:,2))
   rmse_opt_v(i)=sqrt(sum(lrn.sim.indata.data(i,2) - opt.sim.indata.data(i,2))^2/i);
   rmse_naive_v(i)=sqrt(sum(lrn.sim.indata.data(i,2) - naive.sim.indata.data(i,2))^2/i);
end
plot(rmse_opt_v, '.g'); hold on;
plot(rmse_naive_v, '.b'); hold on; box off;
legend(sprintf('Optimizer based decoder recovery - RMSE : %f', rmse_opt),sprintf('Naive decoder recovery - RMSE : %f', rmse_naive));