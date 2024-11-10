% Authors of this file: Mahmoud S. Shaqfa @ 2024

function D_mat = spheroidal_harmonic_basis(max_n, thetas, phis, spheroid_type, varargin)
% Imeplementation: Mahmoud Shaqfa.

if nargin == 4
    printed = true;
else
    printed = varargin{1};
end

if spheroid_type == "oblate"
    xi = sin(thetas);
elseif spheroid_type == "prolate"
    xi = cos(thetas);
end

D_mat = zeros(size(thetas, 1), (max_n+1)^2);

for n = 0:max_n
    Legendre_table = legendre(n, xi)';
%   For positive orders
    for m = 0:n
        % Neumann BCs
        Normalization = sqrt((2*n +1) * factorial(n-m) / (4 * pi * factorial(n+m)));
        D_mat(:, n^2 + n + m + 1) = Normalization .* Legendre_table(:, m+1) .* exp(1i .* m.* phis);
    end
%   For redundunt negative orders
    for m = -n:-1
        D_mat(:, n^2 + n + m + 1) = conj(D_mat(:, n^2 + n + abs(m) + 1)) .* (-1) .^ abs(m);
    end
    
    if printed
        fprintf("\nCompleted: %3.2f %% for n = %3.0f", ((n+1)^2/(max_n+1)^2*100), n)
    end
end