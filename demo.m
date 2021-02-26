clear all; close all;
advanced_throughput = load('./color_advanced/advanced_throughput.mat');
advanced_throughput_var = load('./color_advanced/advanced_throughput_var.mat');

contrast_throughput = load('./color_contrast/contrast_throughput.mat');
contrast_throughput_var = load('./color_contrast/contrast_throughput_var.mat');

figure
plot(5:35, advanced_throughput.throughput_all_rec, '-o')
hold on
plot(5:35, contrast_throughput.throughput_all_rec, '-*')
legend('advanced','contrast','location','SouthEast')


figure
plot(5:35, advanced_throughput_var.throughput_var_rec)
hold on
plot(5:35, contrast_throughput_var.throughput_var_rec)
legend('advanced','contrast')