function [exa, Y2] = hh1_H2(a, b, n, epsilon)
% hh1_H2: Soluciona la ecuación de calor 1D homogénea con condiciones de frontera homogéneas
% mediante una red neuronal basada en series de Fourier.
%
% Entradas:
%   a, b     - extremos del dominio espacial
%   n        - número de nodos espaciales
%   epsilon  - tolerancia para el criterio de parada
%
% Salidas:
%   exa      - solución aproximada (neuronal)
%   Y2       - solución exacta (Fourier)

% ----------------------
% PARTE 1: Configuración
% ----------------------
tic
xs = linspace(a, b, n);    % malla espacial
dt = 1/50;                 % paso temporal
t = dt:dt:1;               % vector de tiempo
nu = 0.5;                  % difusividad térmica
tf = length(t);            % número de pasos temporales
Y = zeros(tf, n);          % solución exacta (intermedia)
Y2 = zeros(n, tf);         % solución exacta (transpuesta para gráficos)

% Funciones fuente
f = @(x) x.*(1 - x);       % fuente inicial
w = @(x) 0;                % condiciones de frontera homogéneas

% ----------------------
% PARTE 2: Solución exacta por series de Fourier
% ----------------------
[~, cn] = Cn_nh(f, w, 1, n);        % coeficientes de Fourier
[X, T] = meshgrid(xs, t);          % rejilla espacio-tiempo

% Construcción de la solución exacta
for i = 1:n
    Y = Y + cn(i) * (sin(i*pi*X) + cos(i*pi*X)) .* exp(-nu*(pi^2)*i^2*T);
end
Y2 = Y';    % transpuesta para coincidir con formato esperado

% ----------------------
% PARTE 3: Aproximación neuronal
% ----------------------
alpha_r = 1e-1;            % tasa de aprendizaje
a_i = zeros(n,1);          % pesos iniciales (neuronales)

for k = 0:1e7
    a_g = zeros(n,1);      % gradiente acumulado

    for j = 1:tf
        yn = 0;

        % Propagación hacia adelante (forward)
        for h = 1:n
            decay = exp(-h^2 * nu * pi^2 * t(j));
            basis = cos(h*pi*X) + sin(h*pi*X);
            yn = yn + tanh(a_i(h) * decay * basis);
        end
        yn = yn';

        % Cálculo del gradiente
        for l = 1:n
            decay = exp(-l^2 * nu * pi^2 * t(j));
            basis_l = cos(l*pi*xs) + sin(l*pi*xs);
            sech2 = (sech(a_i(l)*decay*basis_l)).^2;

            error_term = Y2(:, j) - yn(:, j);
            a_g(l) = a_g(l) - (2/(tf * n)) * sum(error_term .* sech2 .* decay .* basis_l');
        end
    end

    % Actualización de pesos (descenso de gradiente)
    a_i = a_i - alpha_r * a_g;
    disp(alpha_r * max(abs(a_g)));

    % Criterio de convergencia
    if alpha_r * max(abs(a_g)) < epsilon
        toc
        break
    end
end

% ----------------------
% PARTE 4: Reconstrucción de la solución aproximada
% ----------------------
[T, X] = meshgrid(t, xs);
exa = zeros(n, tf);
for i = 1:n
    exa = exa + a_i(i) * (sin(i*pi*X) + cos(i*pi*X)) .* exp(-nu*(pi^2)*i^2*T);
end

% ----------------------
% PARTE 5: Visualización
% ----------------------
figure(1)
mesh(exa)
colormap cool;
title("Numerical solution (Neural approximation)")
view(3)

figure(2)
mesh(Y2)
colormap cool;
title("Exact solution (Fourier)")

figure(3)
mesh(Y2 - exa)
colormap cool;
title("Error (Exact - Approximation)")
view(3)

end
