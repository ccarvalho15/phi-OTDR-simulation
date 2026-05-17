clear; clc; close all

%% =========================================================================
%% SIMULADOR NUMÉRICO DE REFLECTOMETRIA ÓTICA (OTDR / OFDR)
%% =========================================================================
% Este modelo simula a resposta de uma fibra ótica baseada no backscattering
% de Rayleigh. Utiliza perturbações diretas do índice de refração para modelar
% eventos (ex: variação de temperatura/strain) e extrai o desvio espectral
% através da correlação cruzada do speckle ótico.

%% ---------------- 1. PARÂMETROS FÍSICOS E HARDWARE ----------------
c = 3e8;                    % Velocidade da luz no vácuo (m/s)
lambda0 = 1550e-9;          % Comprimento de onda central (m) 
nu0 = c/lambda0;            % Frequência ótica central (~193.5 THz)
n_ave = 1.456;              % Índice de refração médio (Fibra Monomodo - SMF)

L = 1000;                   % Comprimento total da fibra simulada (m)
dz = 0.05;                  % Intervalo de amostragem espacial (m)
z = 0:dz:L-dz;              % Vetor de distância
Nz = length(z);             % Número de pontos espaciais

% --- PARÂMETROS DO PULSO E TRADE-OFF ESPACIAL ---
pulse_width = 10e-9;        % Largura do pulso ótico (ex: 10 ns)

% EQUAÇÃO DA RESOLUÇÃO ESPACIAL: d = (c * pulse_width) / (2 * n_ave)
% O fator '2' deve-se ao percurso de ida e volta do sinal (backscattering).
d = c*pulse_width/(2*n_ave); 

% TRADE-OFF (dz vs. d): 
% Para o backscattering de Rayleigh funcionar, o pulso tem de englobar pontos
% suficientes para gerar um padrão estatisticamente aleatório (speckle/footprint).
% A variável M define o número de scattering virtuais dentro do pulso.
% REGRA: Se dz for muito grande (ex: M < 20), a estatística colapsa,
% o que gera uma onda periódica simples e desformatando a correlação cruzada.
M = round(d/dz);            

fprintf('--- Inicialização do Simulador ---\n');
fprintf('Resolução Física (d): %.2f m | Scattering por pulso (M): %d\n', d, M);

%% ------------ 2. PERFIL DA FIBRA (CENTROS DE RAYLEIGH) -----------
sigma_n = 2e-6;                
delta_n = sigma_n*randn(1,Nz); % Ruído gaussiano simulando as inomogeneidades
n = n_ave + delta_n;           

% EQUAÇÃO DE FRESNEL (Incidência Normal): r = (n1 - n2) / (n1 + n2)
% Calcula a fração de luz refletida em cada fronteira 'dz' da fibra.
r = zeros(1,Nz);
for i = 1:Nz-1
    r(i) = (n(i)-n(i+1))/(n(i)+n(i+1));
end

%% ------------- 3. INJEÇÃO DE EVENTOS E PERTURBAÇÕES ------
delta_n_pert = zeros(1, Nz);

num_events = 5;             
pert_length = 5;            % Comprimento de cada evento (m)
spacing = 150;              % Espaçamento entre eles (m)

% --- TRADE-OFF DE MAGNITUDE (delta_n vs. freq_range) ---
% EQUAÇÃO DO DESVIO: delta_nu / nu0 = - delta_n / n_ave
% O desvio da frequência medida (delta_nu) é proporcional à variação do índice.
% REGRA: Um delta_n de 1e-4 gera desvios de ~13 GHz. Se a nossa janela de 
% observação for de apenas 1 GHz, o pico é atirado para fora do gráfico e 
% o algoritmo falha. Utilizam-se valores na ordem de 10^-7 (dezenas de MHz) 
% para simular efeitos térmicos ou tensões mecânicas mais realistas.
min_mag = 1e-7; 
max_mag = 5e-7;

% O espetro de speckle gera um pico de correlação com uma largura da base 
% (Correlation Bandwidth) a depender do tamanho do pulso na fibra (d), 
% sendo dada por: Largura (FWHM) = c / (2 * n_ave * d)
%
% CÁLCULO PARA 10 ns (d ~ 1 metro):
%    Largura = 3e8 / (2 * 1.456 * 1.03) ~= 100 MHz
%
% Isto significa que o pico de correlação é uma "montanha" com 100 MHz 
% de largura na base. Se um evento gerar um desvio de apenas 20 ou 30 MHz, 
% esse desvio ocorre *dentro* da montanha. O algoritmo 1D encontra a ponta, 
% mas no 3D (surf) a montanha parece apenas ligeiramente deformada.
%
% Novos valores
% - Opção A (Micro-eventos): Se delta_n = 1e-7 -> Desvio de ~19 MHz (Fica dentro da montanha).
% - Opção B (Eventos Destacados): Para o pico sair fora da montanha no 3D,
%   precisamos de desvios > 100 MHz. Isso exige um delta_n > 7.5e-7.
%
% Aqui aplicamos magnitudes ligeiramente mais altas para garantir que os 
% picos saltam para fora da base dos 100 MHz, ficando visualmente muito 
% evidentes no mapa 3D (Desvios entre ~95 MHz e ~190 MHz):

% min_mag = 5e-7; 
% max_mag = 1e-6;


% Cálculo de segurança (previne erros se a geometria exceder a fibra)
total_length = (num_events * pert_length) + ((num_events - 1) * spacing);
if total_length > L
    error('Erro: A configuração dos eventos excede o comprimento da fibra.');
end

random_start = (L - total_length) * rand();
first_event = random_start;
last_event = random_start + total_length;

for i = 1:num_events
    start_idx = max(1, round(random_start/dz));
    end_idx = min(Nz, start_idx + round(pert_length/dz));
    
    % Atribui uma magnitude aleatória dentro dos limites realistas
    random_mag = (min_mag + (max_mag - min_mag) * rand()) * sign(rand - 0.5);
    delta_n_pert(start_idx:end_idx) = random_mag;
    
    random_start = random_start + pert_length + spacing;
end

n_pert = n + delta_n_pert;

%% ------------ 4. MATRIZ DE AQUISIÇÃO E PROPAGAÇÃO ------------
% Simula o varrimento do laser 
freq_range = 1000e6;        % Banda total (1 GHz)
delta_f = 5e6;              % Passo de resolução (5 MHz)
f = nu0 + (-freq_range/2 : delta_f : freq_range/2);  
Nf = length(f);

E_ref = zeros(Nf, Nz-M+1);
E_sig = zeros(Nf, Nz-M+1);
window = ones(1, M);        % O pulso ótico atua como uma janela retangular

fprintf('\nA processar a propagação ótica (Motor Vetorizado)...\n');
%% ------------ FREQUENCY SCANNING LOOP (O MOTOR FÍSICO) ------------
for f_idx = 1:Nf
    current_nu = f(f_idx);
    
    % EQUAÇÃO DA FASE: phi = integral( 2 * pi * n * nu / c ) dz
    % A fase da luz muda consoante a frequência do laser e o índice da fibra.
    phi_ref = cumsum(2 * pi * n * current_nu / c) * dz;
    phi_sig = cumsum(2 * pi * n_pert * current_nu / c) * dz;
    
    % EQUAÇÃO DO CAMPO ELÉTRICO: E(z) = Soma( r(z) * exp(j * 2 * phi) )
    % PORQUÊ CONVOLUÇÃO ('conv')? 
    % Fisicamente, o fotodetector lê a soma fasorial simultânea de todos os 
    % refletores (M) dentro da largura do pulso. Matematicamente, isto é uma 
    % janela deslizante ao longo da fibra.
    % O uso da convolução substitui um loop 'for' espacial encadeado e mais lento
    % por uma operação de matriz nativa. O cálculo torna-se dezenas de vezes 
    % mais rápido produzindo um resultado matemático exatamente igual.
    E_ref(f_idx, :) = conv(r .* exp(1j * 2 * phi_ref), window, 'valid');
    E_sig(f_idx, :) = conv(r .* exp(1j * 2 * phi_sig), window, 'valid');
end

%% ------------ 5. CORRELAÇÃO CRUZADA E DETEÇÃO ------------
fprintf('A extrair desvios de frequência (Cross-Correlation)...\n');

lags_freq = (-(Nf-1):(Nf-1)) * delta_f; 
corr_map = zeros(length(lags_freq), Nz-M+1); 
freq_shift = zeros(1, Nz-M+1); 

for k = 1:Nz-M+1
    % O fotodetetor lê a intensidade (Potência), que é o quadrado do campo complexo
    ref = abs(E_ref(:,k)).^2;
    sig = abs(E_sig(:,k)).^2;

    % Subtração da Média (mean):
    % A intensidade ótica tem um forte componente DC (é sempre positiva). 
    % Subtrair a média limpa o sinal, isolando apenas a "rugosidade" AC (o speckle).
    % Sem isto, a correlação foca-se na intensidade base e não na assinatura de Rayleigh.
    [cv, lags] = xcorr(sig - mean(sig), ref - mean(ref), 'coeff');
    corr_map(:, k) = cv; 

    % O deslocamento do pico máximo indica a variação exata da perturbação
    [~, max_idx] = max(cv);
    freq_shift(k) = lags(max_idx) * delta_f; 
end

%% ------------ 6. VISUALIZAÇÃO DE RESULTADOS ------------
% FIGURA 1: Deteção 1D Simples
figure(1)
plot(z(1:Nz-M+1), freq_shift / 1e6, 'LineWidth', 1.5)
grid on; ylabel('Frequency Shift (MHz)'); xlabel('Distance (m)');
title('1D Sensor Response along the Fiber');
ylim([-250 250]); % Ajustado para os novos limites de magnitude

% FIGURA 2: Mapa 2D com Vista de Topo (Análise Fina / Contorno)
% Como visto em acima, o pico natural tem ~100 MHz de largura. Se testarmos
% micro-tensões (< 100 MHz), elas ficam visualmente engolidas na inclinação 
% da montanha tridimensional. A vista de topo atua como um mapa 
% topográfico: corta a montanha em fatias e permite ver o interior do s
% peckle ótico a desviar-se limpidamente, por menor que seja o desvio.

figure(2)
[Z_mesh, F_mesh] = meshgrid(z(1:Nz-M+1), lags_freq / 1e6);
contourf(Z_mesh, F_mesh, corr_map, 20, 'LineColor', 'none'); 
colormap('jet'); colorbar;
hold on;
% A linha branca traça o pico detetado pelo algoritmo matemático por cima do radar ótico
plot(z(1:Nz-M+1), freq_shift / 1e6, 'w', 'LineWidth', 1.5); 
hold off;
xlabel('Distance (m)'); ylabel('Frequency Lag (MHz)');
title('2D Correlation Map (Top View with Peak Trace)');
xlim([max(0, first_event - 50) min(L, last_event + 50)]); 
ylim([-250 250]);

% FIGURA 3: Visão Integrada (1D e 3D)
figure(3)
subplot(2,1,1)
plot(z(1:Nz-M+1), freq_shift / 1e6, 'LineWidth', 1.5)
grid on; ylabel('Frequency Shift (MHz)'); xlabel('Distance (m)');
title('Total Disturbance Profile (Frequency Shift)');
ylim([-250 250]);

subplot(2,1,2)
surf(Z_mesh, F_mesh, corr_map, 'EdgeColor', 'none')
view(35, 45); colormap('jet'); colorbar;
xlabel('Distance (m)'); ylabel('Frequency Lag (MHz)'); zlabel('Correlation');
title('3D Correlation Surface (Full Fiber)');
rotate3d on;
xlim([0 L]); 
ylim([-250 250]); 
zlim([-0.5 1]);

fprintf('--- Simulação concluída com sucesso! ---\n');