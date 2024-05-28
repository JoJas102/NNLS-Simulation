function [dNNLS, fNNLS, resultNNLS, dAUC, fAUC] = findpeaksNNLS(matrix, np, reg)
    
    % fit with (reg = 1) and w\ (= 0) regression
    if reg == 1 % sNNLSReg
        s(:,:) =  matrix(3,:,:); 
    else        % sNNLSNoReg
        s(:,:) =  matrix(2,:,:); 
        fprintf('ATTENTION: regularisation of NNLS results is deactivated \n')
    end
    
    DValues =  matrix(1,:,1);
    
    % right array dimensions and #iter for single simu
    if ndims(matrix) < 3 %#ok<ISMAT> 
        iter = 1;
        s = s'; 
    else
        iter = length(matrix(1,1,:));
    end
    
    % Matrix indices at iteration i
    % [1,i] slow/tissue
    % [2,i] inter/tubular
    % [3,i] fast/blood
    
    [dNNLS,fNNLS] = deal(zeros(np+1, iter));
    resultNNLS = zeros(np+1,4,iter);
    [fAUC, dAUC] = deal(zeros(iter,3));
    
    for i=1:iter
        clear maxima d widths proms;
        fp = length(findpeaks(s(:,i),'MinPeakHeight',0.06));                    % found number of peaks
        
        [maxima(:),d(:),widths(:), ~] = findpeaks(s(:,i), 1:length(DValues),... % returns maxima and their x positions as indices plus the FWHM and 
            'WidthReference','halfheight','MinPeakHeight',0.06, 'NPeaks',np+1); % 'SortStr','descend', % prominences in descending peak height
        %findpeaks(s(:,i), 1:length(DValues), 'Annotate','extents',...        % check findpeaks visually
        %'WidthReference','halfheight','MinPeakHeight',0.06,'NPeaks',np+1); 
                
        d = DValues(d) ;                                                        % convert back to log scale values
        f = maxima.*widths/(2*sqrt(2*log(2)))*sqrt(2*pi);                       % calc area under gaussian curve (= vol frac)       
        f = f./sum(f)*100;                                                      % normalise f
        result = [maxima; d; widths; f]';
        
        if fp < np+1                                                            % fill with zeros if < np+1 peaks are found
            nz = np+1 - fp;
            result = [result; zeros(nz,4)];                                     % result(end+1:numel(A))=0
        end
        
        resultNNLS(:,:,i) = result;
        dNNLS(:,i) = result(:,2);
        fNNLS(:,i) = result(:,4);
        
        % AUC calculation
        [dAUC(i,:), fAUC(i,:)] = findAUC(s, dNNLS, fNNLS, DValues, i);

        % Threshold
        dNNLS(fNNLS<3)=0;                                                       % remove any entry with vol frac lower 3%
        fNNLS(fNNLS<3)=0;
        
        if fNNLS(4,i) ~= 0 && fNNLS(3,i) == 0                                   % check for threshhold error/fourth peak and fill up
            fNNLS([3 4],i) = fNNLS([4 3],i);                                    % HARD CODED!
            dNNLS([3 4],i) = dNNLS([4 3],i);
        end
        
        [fNNLS(:,i), sortID] = sort(fNNLS(:,i),'descend');                      % sort array for descending vol fraction (4th entry = smallest peaks)
        dNNLS(:,i) = dNNLS(sortID,i);

        fNNLS(:,i) = fNNLS(:,i)./sum(fNNLS(:,i))*100;                           % normalise again after filtering for too low peaks

        resultNNLS(:,2,i) = dNNLS(:,i);                                                  % update resultNNLS 
        resultNNLS(:,4,i) = fNNLS(:,i);
    end
        
end