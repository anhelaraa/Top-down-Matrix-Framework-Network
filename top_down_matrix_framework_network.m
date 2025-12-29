clc;
clear;
close all

% =========================================================================
% This code solves a NxN linear system Ax = b using Gaussian eliminations,
% with scaled partial pivoting to swap rows when needed for stability,
% checks pivots against a tolerance to detect singular/near-singular
% systems, then finishes with back substitution to compute the solution
% vector x. If the solve succeeds it prints x; otherwise it reports an
% error.      
% =========================================================================

A = [2 -1 0;
    -1 4 -3;
    0 -3 9];
b = [12; -4; 0];

n = 3;
tol = 1e-12;

[x, er] = gauss(A, b, n, tol); % Call gauss function

% Display the results
if er == 0
    fprintf('The solution vector for x is:\n');
    format bank; % Only shows 2 decimal places.
    disp(x);
else
    fprintf('An error occurred. The system is singular or near-singular.\n');
end

% ======================================================
%                      Functions
% ======================================================
function [a, b, er] = eliminate(a, b, s, n, tol, er)

for k = 1:n - 1
    [a, b, s, er] = pivot(a, b, s, n, k, er); % CALL Pivot
    if er == -1
        return;
    end

    if abs(a(k, k) / s(k)) < tol
        er = -1;
        break; % EXITDO
    end

    for i = k + 1:n
        factor = a(i, k) / a(k, k);
        for j = k:n
            a(i, j) = a(i, j) - factor * a(k, j);
        end
        b(i) = b(i) - factor * b(k);

    end
end
end

function [a, b, s, er] = pivot(a, b, s, n, k, er) % SUB Pivot
p = k;
big = abs(a(k, k) / s(k)); % s from gauss

for i = k + 1:n % ii
    dummy = abs(a(i, k) / s(i));
    if dummy > big
        big = dummy;
        p = i;
    end
end

% If a swap is needed, perform the row swap
if p ~= k
    % Swap rows in matrix A
    dummy_row_a = a(p, :);
    a(p, :) = a(k, :);
    a(k, :) = dummy_row_a;

    % Swaps the numbers in vector b
    dummy_b = b(p);
    b(p) = b(k);
    b(k) = dummy_b;

    % Swaps the numbers in scaling vector s
    dummy_s = s(p);
    s(p) = s(k);
    s(k) = dummy_s;
end
end

function x = substitute(a, b, n)
x = zeros(n, 1);

x(n) = b(n) / a(n, n);

for i = n - 1:-1:1
    sum = 0;
    for j = i + 1:n
        sum = sum + a(i, j) * x(j);
    end
    x(i) = (b(i) - sum) / a(i, i);
end
end