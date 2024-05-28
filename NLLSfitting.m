function [dNLLS, fNLLS, resnorm] = NLLSfitting(inputSimu, Dmin, Dmax, dIn, fIn, algorithm)
% NLLSfitting(inputSimu, Dmin, Dmax, dIn, fIn) =  a priori information dNNLS and fNNLS in x0
% NLLSfitting(inputSimu, Dmin, Dmax,) = no a priori information, using standard start value

    numberOfIter = length(inputSimu(1,1,:));
    options.Algorithm = 'trust-region-reflective'; 
    options.Display = 'off';
    %options.StepTolerance = 1e-20;
    %options.FunctionTolerance = 1e-10;
    %options.Diagnostics = 'on';
    
    if nargin < 4 % default start values for triExp [Periquito2021]
        %input = repelem([1.5*1e-3, 10*1e-3, 140*1e-3, 56, 36, 8]',1,numberOfIter);
        input = repelem([1.5*1e-3, 30*1e-3, 100*1e-3, 55, 25, 20]',1,numberOfIter); %bad start values [Periquito2021]&Stabinska
        %input = repelem([1.5*1e-3, 30*1e-3, 500*1e-3, 55, 25, 20]',1,numberOfIter); % test
        dIn = input(1:3,:);
        fIn = input(4:6,:);
    elseif nargin == 6
        options.Algorithm = algorithm;
    end
    
    input = [dIn; fIn];
    b = inputSimu(1,:,1);
    np = nnz(inputSimu(2,:,1));
    resultNLLS = zeros(2*(np+1),numberOfIter);
    resnorm = zeros(1,numberOfIter);

    

    for i = 1:numberOfIter
        clear result;
        printOut = fprintf('Simulating NLLS... %.1f%% \n', i/numberOfIter*100);
        
        fp = nnz(fIn(:,i));                       % number of found compartments (by NNLS)
        lb = [repelem(Dmin,fp) repelem(0,fp-1)];  % set bound constraints based on NNLS d range
        ub = [repelem(Dmax,fp) repelem(100,fp-1)];
        signal(:,:) = inputSimu(5,:,i);
        x0 = nonzeros(input(:,i));                % start values without zeros for fitting
        x0 = x0(1:end-1,:);                       % x0 = [d_1 ... d_N f_1 ... f_N-1]
        
        if fp == 3 
            % Create tri-exponential signal function for fitting with d and f as fitting variable
            triExp = @(x) exp(-kron(b, abs(x(1))))*x(4) + exp(-kron(b, abs(x(2))))*x(5) + exp(-kron(b, abs(x(3))))*(100-(x(4)+x(5))) - signal;
            [result, resnorm(:,i)] = lsqnonlin(triExp, x0, lb, ub, options); % options
            result(6) =(100-(result(4)+result(5)));
            resultNLLS(:,i) = [result(1:3); zeros(np+1-fp,1); result(4:6); zeros(np+1-fp,1)];
            
        elseif fp == 2 
            % Create bi-exponential signal function
            biExp = @(x) exp(-kron(b, abs(x(1))))*x(3) + exp(-kron(b, abs(x(2))))*(100-x(3)) - signal;
            [result, resnorm(:,i)] = lsqnonlin(biExp, x0, lb, ub, options); 
            result(4) =(100-result(3));
            resultNLLS(:,i) = [result(1:2); zeros(np+1-fp,1); result(3:4); zeros(np+1-fp,1)];
            
        elseif fp == 1 
            % Create mono-exponential signal function
            monoExp = @(x) exp(-kron(b, abs(x(1)))) - signal;
            [result, resnorm(:,i)] = lsqnonlin(monoExp, x0, lb, ub, options); 
            result(2) = 1;        
            resultNLLS(:,i) = [result(1); zeros(np+1-fp,1); result(2); zeros(np+1-fp,1)];
            
        elseif fp == 4 
            % Create 4-exponential signal function
            quadExp = @(x) exp(-kron(b, abs(x(1))))*x(5) + exp(-kron(b, abs(x(2))))*x(6) + exp(-kron(b, abs(x(3))))*x(7) + exp(-kron(b, abs(x(4))))*(100-(x(5)+x(6)+x(7))) - signal;
            [result, resnorm(:,i)] = lsqnonlin(quadExp, x0, lb, ub, options); 
            result(8) =(100-(result(5)+result(6)+result(7)));
            resultNLLS(:,i) = result;
        end
        dNLLS(:,i) = abs(resultNLLS(1:np+1,i));
        fNLLS(:,i) = resultNLLS(np+2:end,i);

        [fNLLS(:,i), sortID] = sort(fNLLS(:,i),'descend');                      % sort array for descending vol fraction (4th entry = smallest peaks)
        dNLLS(:,i) = dNLLS(sortID,i);

        for j = 1:printOut; fprintf('\b'); end
    end
    
    %dNLLS = abs(resultNLLS(1:np+1,:));
    %fNLLS = resultNLLS(np+2:end,:);

end