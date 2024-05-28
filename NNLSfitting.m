function [s, s2, mu] = NNLSfitting(DBasis, signal, i)

    printOut = fprintf('Simulating NNLS... %.1f%%', i);
    
    % NNLS fitting w\o Reg (find unknown signal amplitude of components (signal)) minimises norm(A*s-signal) for reference
    s = lsqnonneg(DBasis,signal); 

    % Regularization fitting NNLS with CVNNLS from Bjarnason
    [s2, mu, resid] = CVNNLS(DBasis,signal);
    
    %fprintf('µ = %f \n', mu); % larger mu = more satisfaction of constraints at expanse of increasing misfit (Witthal1989)

    for j = 1:printOut; fprintf('\b'); end
end