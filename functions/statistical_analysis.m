function [Precision, Sensitivity, F1] = statistical_analysis(gt, conf, ...
    all_locs)
    %% ===================================================================
    % 8.2 STATISTICAL PERFORMANCE ANALYSIS (TP, FP, FN)
    % ====================================================================
    %   Performs spatial matching between ground-truth configurations and
    %   experimental sensor peaks to classify detections into a Confusion 
    %   Matrix.
    %       - True Positives (TP): Real perturbations successfully resolved.
    %       - False Negatives (FN): Real perturbations missed or masked by 
    %       noise.
    %       - False Positives (FP): Noise spikes falsely classified as 
    %       events.

    match_tolerance = conf.match_tolerance;

    num_events = length(gt.theor_starts);
    num_detected = length(all_locs);
    % Track which detected peaks match actual ground truth

    detected_matched = false(1, num_detected); 
    TP = 0; % True Positives counter (correctly identified events)
    FN = 0; % False Negatives counter (missed events)
    
    % Evaluate TP and FN against ground truth
    for i = 1:num_events
        % Define spatial acceptance boundaries using the tolerance margin
        limit_start = gt.theor_starts(i) - match_tolerance;
        limit_end = gt.theor_ends(i) + match_tolerance;
        
        % Check if any detected peaks fall within the valid spatial window
        match_idx = find(all_locs >= limit_start & all_locs <= limit_end);
        
        if ~isempty(match_idx)
            theor_center = (gt.theor_starts(i) + gt.theor_ends(i)) / 2;

            % [~, local_best] = min(abs(all_locs(match_idx) - theor_center));
            [max_val, local_best] = max(abs(all_pks(match_idx)));
            best_peak_idx = match_idx(local_best);
            
            if max_val >= threshold
                TP = TP + 1;
                detected_matched(best_peak_idx) = true;
            else
                FN = FN + 1;
            end
        else
            FN = FN + 1;
        end
    end
    
    % Evaluate FP: they're detected peaks that did not map to any real 
    % ground truth event
    FP = sum(~detected_matched);
    
    % Calculate key signal detection metrics: Sensitivity, Precision, and 
    % F1-Score
    Sensitivity = (TP / (TP + FN)) * 100; 
    Precision = 0;
    F1 = 0;
    if (TP + FP) > 0
        Precision = (TP / (TP + FP)) * 100;
    end
    if (Precision + Sensitivity) > 0
        F1 = 2 * (Precision * Sensitivity) / (Precision + Sensitivity);
    end

%% 
% TP (True Positive) = O sensor detetou um pico e existia lá uma perturbação real
% FP (False Positive) = O sensor detetou um pico, mas era apenas ruído (falso alarme)
% FN (False Negative) = Existia uma perturbação real, mas o sensor não a detetou (falha)

% Cálculo de precisão
% A precisão responde à pergunta: "De todos os alarmes que o sensor disparou, 
% quantos eram reais?"
% Ela foca-se na fidelidade do sinal. Se tiveres uma precisão de 100%, 
% significa que nunca tens alarmes falsos, mas não garante que detetaste 
% todos os eventos.
% FÓRMULA :: prec = TP / (TP + FP)
%            prec = TP / max((TP + FP), 1);

% Cálculo da sensibilidade
% A sensibilidade (também chamada de Recall) responde à pergunta: "De todos 
% os eventos que aconteceram na fibra, quantos é que eu consegui apanhar?"
% Ela foca-se na abrangência. Se tiveres um Recall de 100%, significa que 
% não deixaste escapar nenhuma perturbação, mas podes ter disparado muitos 
% alarmes falsos pelo caminho.
% FÓRMULA :: rec = TP / (TP + FN)
%            rec = TP / max((TP + FN), 1);

% Cálculo da F1-Score
% O F1-Score é a média harmónica entre a Precisão e a Sensibilidade. É a 
% métrica mais usada em artigos científicos (como o de Zheyuan Lin que 
% carregaste) porque resume o desempenho num único número.
% 
% Porquê usar a média harmónica e não a média normal?
% Porque a média harmónica penaliza valores extremos. Se a tua Precisão for 
% 1.0 (perfeita) mas o teu Recall for 0.0 (não detetaste nada), a média 
% normal seria 0.5, mas o F1-Score será próximo de 0. O F1-Score só é alto 
% se ambas as métricas forem boas.
% FÓRMULA :: f1 = 2 * (prec * rec) / (prec + rec)
%            F1 = (2 * prec * rec) / (prec + rec);


% EXEMPLO PRATICO DE ANALISE
% Simulação para 5 eventos reais
% - O código deteta 6 picos
% - Ao comparar as distâncias, vês que 4 picos batem certo com os eventos 
% reais (TP = 4).
% - 2 picos eram apenas ruído que passou a threshold (FP = 2)
% - 1 evento real foi ignorado pelo sensor (FN = 1)
% RESULTADOS:
% - prec = 66% dos alarmes são verdadeiros
% - rec = 80% eventos reais
% - f1 = 72%

% OBJETIVO: encontrar o ponto onde este valor chega o mais perto possível 
% de 1.0.