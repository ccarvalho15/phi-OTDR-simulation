% 0. Define custom normpdf (bypasses Statistics Toolbox requirement)
normpdf = @(x, mu, sigma) (1 / (sigma * sqrt(2*pi))) * exp(-((x - mu).^2) / (2 * sigma^2));

% 1. Define distribution parameters
mu = 0;      % Mean
sigma = 1;   % Standard deviation

% 2. X vector and PDF calculation
x = linspace(mu - 4*sigma, mu + 4*sigma, 1000);
y = normpdf(x, mu, sigma);

figure('Color', 'w');
hold on;

% 3. Fill regions
x3 = linspace(mu - 3*sigma, mu + 3*sigma, 500);
fill([x3, fliplr(x3)], [normpdf(x3, mu, sigma), zeros(1, length(x3))], ...
     [0.90 0.93 1.00], 'EdgeColor', 'none', 'DisplayName', '\pm 3\sigma (99.73%)');

x2 = linspace(mu - 2*sigma, mu + 2*sigma, 500);
fill([x2, fliplr(x2)], [normpdf(x2, mu, sigma), zeros(1, length(x2))], ...
     [0.70 0.80 0.98], 'EdgeColor', 'none', 'DisplayName', '\pm 2\sigma (95.45%)');

x1 = linspace(mu - 1*sigma, mu + 1*sigma, 500);
fill([x1, fliplr(x1)], [normpdf(x1, mu, sigma), zeros(1, length(x1))], ...
     [0.45 0.62 0.95], 'EdgeColor', 'none', 'DisplayName', '\pm \sigma (68.27%)');

% 4. Plot main curve
plot(x, y, 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% 5. Add vertical reference lines
xline(mu, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
for s = -3:3
    if s ~= 0
        xline(mu + s*sigma, 'k:', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
    end
end

set(gca, 'FontSize', 16); % Define o tamanho dos números de X e Y

% 6. Figure styling
title('Empirical Rule', 'FontSize', 20,'FontWeight', 'bold');
xlabel('Standard Deviations from the Mean (\sigma)', 'FontSize', 18);
ylabel('Probability Density','FontSize', 18);
legend('Location', 'northeast','FontSize',16);
grid on;
box on;
hold off;