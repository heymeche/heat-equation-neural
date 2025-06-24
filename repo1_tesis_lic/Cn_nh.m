function [s, cn] = Cn_nh(f, w, L, n)
% Cn_nh: Calcula la serie de Fourier con condiciones de frontera homogéneas
% para la ecuación de calor en el intervalo [0, L]
%
% Entradas:
%   f - función fuente inicial (u(x,0))
%   w - función de frontera (no usada aquí, solo compatibilidad)
%   L - longitud del dominio (generalmente L = 1)
%   n - número de términos en la serie
%
% Salidas:
%   s  - aproximación de la función f(x) usando n términos
%   cn - coeficientes cn para cada término seno

cn = zeros(n,1);
s = @(x) 0;

for i = 1:n
    % Definición del término ortogonal base
    phi_i = @(x) sin(i * pi * x / L);

    % Integral para calcular el coeficiente cn(i)
    % Fórmula: cn(i) = 2/L * ∫ f(x) * sin(iπx/L) dx de 0 a L
    cn(i) = (2 / L) * integral(@(x) f(x) .* phi_i(x), 0, L);

    % Suma parcial de la serie
    s = @(x) s(x) + cn(i) * phi_i(x);
end

end
