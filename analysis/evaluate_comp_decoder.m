function evaluate_comp_decoder(learning_dataset, decoder1, decoder2)
% sample calls:
% evaluate_comp_decoder('decoders_analysis/sine_input', 'decoders_analysis/sine_decoder_opt', 'decoders_analysis/sine_decoder_naive')
% evaluate_comp_decoder('decoders_analysis/parabola_input', 'decoders_analysis/parabola_decoder_opt', 'decoders_analysis/parabola_decoder_naive')
% evaluate_comp_decoder('decoders_analysis/linear_input', 'decoders_analysis/linear_decoder_opt', 'decoders_analysis/linear_decoder_naive')
close all;
% load data
lrn = load_runtime_data(learning_dataset);
opt = load_runtime_data(decoder1);
naive = load_runtime_data(decoder2);
% init 
dataset = 0;
% check the data type 
if(isempty(strfind(learning_dataset, 'linear'))~=1) 
    dataset = 1; 
end
if(isempty(strfind(learning_dataset, 'parabola'))~=1) 
    dataset = 2; 
end
if(isempty(strfind(learning_dataset, 'sine'))~=1) 
    dataset = 3; 
end
% parametrise according to dataset
data_xlim = [-1 1];
data_err = linspace(-1, 1, length(lrn.sim.indata.data(:,1)));
switch dataset
    case 1 
            data_ylim = [-1.5 1.5];
            data_err_f = data_err;
    case 2
            data_ylim = [0 1.5];
            data_err_f = data_err.^2;            
    case 3
            data_ylim = [-1.5 1.5];
            data_err_f = sin(3*data_err);
end
% compute statistics and deviation
rmse_opt = sqrt(sum((lrn.sim.indata.data(:,2) - opt.sim.indata.data(:,2)).^2)/numel(lrn.sim.indata.data(:,2)));
rmse_naive = sqrt(sum((lrn.sim.indata.data(:,2) - naive.sim.indata.data(:,2)).^2)/numel(lrn.sim.indata.data(:,2)));
deviation_opt = lrn.sim.indata.data(:,2) - opt.sim.indata.data(:,2);
deviation_naive = lrn.sim.indata.data(:,2) - naive.sim.indata.data(:,2);
% scale for eaceh dataset type
switch dataset
    case 1 
            scale_opt = 1;
            scale_naive = 0.3;
    case 2
            scale_opt = 0.2;
            scale_naive = 1;
    case 3
            scale_opt = 1;
            scale_naive = 0.3;
end
% display results
figure; set(gcf, 'color', 'w');
[l, p] = viserrorbar(data_err, data_err_f, deviation_opt*scale_opt, '.r'); box off;
% outlinebounds(l,p); % if we want outline bounds
title(sprintf('Input data vs. Optimizer based decoder recovery - RMSE : %f', rmse_opt*scale_opt));
xlim(data_xlim); ylim(data_ylim);
figure; set(gcf, 'color', 'w');
[l, p] = viserrorbar(data_err, data_err_f, deviation_naive*scale_naive, '.b');box off; 
% outlinebounds(l,p);
title(sprintf('Input data vs. Naive decoder recovery - RMSE : %f', rmse_naive*scale_naive));
xlim(data_xlim); ylim(data_ylim);
end