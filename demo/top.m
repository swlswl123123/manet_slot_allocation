clear all; close all;
% rng(0)
N = 16;
M = 50;
R = 1000;
t_last = 10;

[advanced_trans, advanced_delay] = data_trans_advanced(N, M, R, t_last);
[color_trans, color_delay] = data_trans_color(N, M, R, t_last);
[random_trans, random_delay] = data_trans_random(N, M, R, t_last);

figure
plot(advanced_trans, '-o')
hold on;
plot(color_trans, '-*')
hold on;
plot(random_trans, '--')

legend('advanced', 'colored', 'random')

length(find(advanced_delay <= 1))
length(find(color_delay <= 1))
length(find(random_delay <= 1))
