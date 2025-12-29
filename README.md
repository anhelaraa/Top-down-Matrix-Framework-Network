# Top-Down Matrix Framework Network

This code solves a NxN linear system Ax = b using Gaussian elimination,
% with scaled partial pivoting to swap rows when needed for stability,
% checks pivots against a tolerance to detect singular/near-singular
% systems, then finishes with back substitution to compute the solution
% vector x. If the solve succeeds, it prints x; otherwise, it reports an
% error.      
