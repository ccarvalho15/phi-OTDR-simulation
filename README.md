# 📡 Phi-OTDR Multi-Perturbation Simulator

Este repositório contém um framework avançado em **MATLAB** para a simulação de sistemas de **Sensoriamento Distribuído em Fibra Óptica** ($\phi$-OTDR). O diferencial deste simulador é a capacidade de modelar e detectar simultaneamente perturbações de **Temperatura** e **Deformação Mecânica (_Strain_)**.

## Destaques do Projeto

O sistema utiliza a técnica de varredura de frequência e correlação cruzada para extrair o deslocamento espectral causado por mudanças físicas na fibra.

* **Modelagem Física Fiel:** Inclui o cálculo de coeficientes termo-ópticos e de expansão térmica.
* **Detecção de Deformação (Strain):** Simula o alongamento/compressão da fibra e seu impacto no índice de refração.
* **Análise Espectral:** Implementa a reconstrução do sinal através da correlação de espectros Rayleigh com média zero (*zero-mean*).
* **Visualização Multidimensional:** Mapas 3D de correlação que facilitam a localização espacial dos eventos.

## Organização dos Módulos

### 1. `static_n.m`
O script principal que integra as duas frentes de perturbação:
* **Perturbação Térmica:** Baseada no coeficiente termo-óptico ($\gamma$) e expansão térmica ($\eta$).
* **Perturbação de Deformação:** Modela a alteração do índice de refração e da fase devido ao estresse mecânico longitudinal na fibra (Strain).

### 2. `refractive_index.m`
Focado na análise estatística do índice de refração sob influência de eventos térmicos aleatórios, permitindo configurar múltiplos eventos ao longo de quilómetros de fibra e visualizar os resultados em diferentes perspectivas 3D.

### 3. `refractive_n_index.m`
Script de validação rápida para perturbações diretas no índice de refração ($\Delta n$), ideal para testes de sensibilidade do algoritmo de detecção sem a complexidade das variáveis ambientais.

## Princípios Físicos

A mudança na frequência óptica ($\Delta \nu$) detectada é proporcional às variações externas conforme a relação fundamental:

$$\frac{\Delta \nu}{\nu} = -\left( K_T \Delta T + K_\epsilon \epsilon \right)$$

Onde:
* $K_T$ é a sensibilidade térmica.
* $K_\epsilon$ é o coeficiente de sensibilidade à deformação (*strain*).
* $\epsilon$ é a deformação aplicada.

## Requisitos e Uso

1.  **Software:** MATLAB R2020a ou superior.
2.  **Toolboxes:** *Signal Processing Toolbox*.
3.  **Execução:**
    * Abra o MATLAB na pasta do projeto.
    * Execute o arquivo desejado (ex: `run('static_n.m')`).
    * Siga as instruções no console para definir o comprimento da fibra e o número de eventos.

## Visualização
Os scripts geram mapas de correlação cruzada onde o pico indica a posição exata da perturbação e a magnitude do deslocamento de frequência (MHz), permitindo o monitoramento em tempo real da integridade estrutural da fibra.
