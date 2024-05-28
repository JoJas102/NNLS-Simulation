function [dAUC, fAUC] = findAUC(s, dNNLS, fNNLS, DValues, i)

    % Find intervals
    [~, index(1)] = min(abs(DValues-0.002));
    [~, index(2)] = min(abs(DValues-0.05));

    % AUC for signal s
    fAUC(1) = trapz(1:index(1),s(1:index(1),i));                  % to compare NNLS results to AUC (>4 peaks)
    fAUC(2) = trapz(index(1):index(2),s(index(1):index(2),i));
    fAUC(3) = trapz(index(2):length(s(:,i)),s(index(2):end,i));
    fAUC(:) = fAUC(:)./sum(fAUC(:))*100;                          % normalise fAUC
    dAUC(1) = dNNLS(1,i);
    dAUC(2) = dNNLS(2,i);
    dAUC(3) = (fNNLS(3,i)*dNNLS(3,i)+fNNLS(4,i)*dNNLS(4,i))/(fNNLS(3,i)+fNNLS(4,i));
    dAUC(isnan(dAUC)) = 0;                                        % 2 peaks save
    %dAUC(3) = dNNLS(3); % testing for 3 peaks
end