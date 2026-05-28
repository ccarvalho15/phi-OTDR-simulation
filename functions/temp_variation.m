function delta_T = temp_variation (delta_n, gamma, eta, n_ave)
    % De acordo com o artigo, a variação equivalente do índice de refração
    % induzida termicamente (dada por delta_n_T) é calculada como
    % delta_n_T = (gamma + n_ave * eta) * delta_T
    % 
    % Para encontrar a variação da temperatura:
    % delta_T = (delta_n_T) / (gamma + n_ave * eta)
    % onde os parâmetros (considerando os valores de referência das fontes
    % para uma fibra padrão a 1550 nm e 300 K) são
    % - Coeficiente termo-ótico (gamma), representa a mudança do índice de
    % refração com a temperatura (~9.1 x 10^-6 K^-1)
    % - Coeficiente de expansão térmica (eta), representa a mudança física
    % no comprimento da fibra (~0.5 x 10^-6 K^-1)
    % - Índice de refração médio (n_ave)

    den = gamma + (n_ave .* eta);
    delta_T = delta_n ./ den;
end