%% بخش اول: تولید سیگنال (این بخش باید حتماً قبل از حلقه باشد)
N = 1000;
rho = 0.05;
W = (rand(1,N) < rho) .* randn(1,N); % تولید نویز خلوت
alpha = 0.8;
% تولید سیگنال خروجی Y
Y = filter(1, [1, -alpha], W); 

%% بخش دوم: تخمین بعد با حلقه For
eps_values = [0.1, 0.05, 0.02, 0.01]; 
results_d = zeros(size(eps_values)); 
results_dW = zeros(size(eps_values)); 
alpha_renyi = 2; 

for i = 1:length(eps_values)
    current_eps = eps_values(i);
    
    % تخمین برای خروجی Y (حالا Y تعریف شده است)
    edges_Y = min(Y):current_eps:max(Y);
    counts_Y = histcounts(Y, edges_Y);
    p_Y = counts_Y / sum(counts_Y);
    p_Y = p_Y(p_Y > 0);
    H_Y = (1 / (1 - alpha_renyi)) * log2(sum(p_Y.^alpha_renyi));
    results_d(i) = H_Y / log2(1/current_eps);
    
    % تخمین برای ورودی W
    edges_W = min(W):current_eps:max(W);
    counts_W = histcounts(W, edges_W);
    p_W = counts_W / sum(counts_W);
    p_W = p_W(p_W > 0);
    H_W = (1 / (1 - alpha_renyi)) * log2(sum(p_W.^alpha_renyi));
    results_dW(i) = H_W / log2(1/current_eps);
end

%% بخش سوم: نمایش نتایج
T = table(eps_values', results_dW', results_d', 'VariableNames', {'Epsilon', 'Dim_Input_W', 'Dim_Output_Y'});
disp('--- Professional Dimension Estimation Table ---');
disp(T);

% رسم نمودار
figure('Color', 'w');
plot(eps_values, results_d, '-ro', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on;
plot(eps_values, results_dW, '-bo', 'LineWidth', 2, 'MarkerFaceColor', 'b');
set(gca, 'XDir', 'reverse'); 
grid on;
legend('Output Signal (Y)', 'Input Noise (W)');
xlabel('Epsilon (Resolution)');
ylabel('Information Dimension (d)');
title('Stability of d Estimation as Epsilon goes to 0');
%% بخش چهارم: رسم مجدد سیگنال‌های ورودی و خروجی
figure('Color', 'w', 'Name', 'Time Domain Signals');

% رسم نویز ورودی خلوت
subplot(2,1,1);
stem(W(1:200), 'Color', [0.4 0.4 0.4], 'MarkerFaceColor', 'b');
title('Input Sparse Noise (W)');
ylabel('Amplitude');
grid on;

% رسم سیگنال خروجی ARMA
subplot(2,1,2);
plot(Y(1:200), 'LineWidth', 1.5, 'Color', 'r');
title('Output ARMA Signal (Y)');
ylabel('Amplitude');
xlabel('Sample Number');
grid on;
% تنظیم محدوده محورها برای مقایسه بهتر
subplot(2,1,1); ylim([-2.5 0.5]);
subplot(2,1,2); ylim([-2.5 0.5]);