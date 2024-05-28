function [d,f,fdInput,DBasis,rawSignal] = createSignal(d, f ,b , Dmin, Dmax, m)

    n = length(b);           %np = length(dRange); % number of input peaks

    % Generate basis values           
    DValues = logspace(log10(Dmin), log10(Dmax), m); 
    DBasis = exp(-kron( b', DValues));
    
    % Real d-values (exact bins)
    [dMin(1), index(1)] = min(abs(DValues-d(1)));
    [dMin(2), index(2)] = min(abs(DValues-d(2)));
    [dMin(3), index(3)] = min(abs(DValues-d(3)));
    d = [DValues(index(1)) DValues(index(2)) DValues(index(3))];
    fdInput = [f zeros(1,n-length(f)); d zeros(1,n-length(d))];
    
    % Calculate offset to reference values for evaluation when randomizing gT
    %fd_offset(:,:) = [d - d_ref 0;
    %             f - f_ref 0];

    % Define total diffusion signal decay according to sum of f_i*e^(-b*D_i)
    A = exp(-kron(b', d)); % constraint matrix containing exp decay funxtions
    rawSignal = A*f';      % mean signal of ROI/voxel without noise
end