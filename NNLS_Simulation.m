%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          NNLS Simulation
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Program to simulate influence of different parameters on NNLS fitting of
% diffusion data and compare multi-exponential fitting methods. Evaluating 
% the reliability of finding the total number of components contributing to
% the corresponding multi-exponantial signal and analysing the results by
% calculating corresponding diffusion parameters, comparing those to the
% ground truth.
% This fitting routine uses the regularized NNLS algorithm with cross 
% validation from Thorarin Bjarnason for comparison.
%
% To carry out this simulation the following functions need to be located
% in the MATLAB path:
% InitVar.m
% DiffParamRandomizer.m
% createNoise.m
% NNLSfitting.m 
% CVNNLS.m
% fastnnls.m
% findpeaksNNLS.m
% findAUC.m
% NLLSfitting.m
% evaluation.m
% plotSimu.m 
% write3DMatrixToTxt.m 
%
% Simulation results and diffusion components estimates will be saved in 
% the result subfolder.
%
% For any further comments or explaining descriptions see annotations
% inside the functions code files.
%
% --
% Jonas Jasse 
% Last modified 16.11.2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial variables mandatory (InitVar.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% - numberOfIter = number of simulations/iteration steps, differentiation 
%                  for single- or multi-plotting of results
% - snrIn        = Input SNR for artificial DWI signal rawSignal()
% - b[]          = array of b-values
% - Dmin/Dmax    = defining the D-range for diffusion constants/basis
%                  values DValues[] and DBasis[]
% - m            = number of diffusion coefficients/bins between Dmin and Dmax
%
% eventually:
% - f[] and d[]  = arrays with volume fractions [tissue tubules blood]/ 
%                  [slow, inter and fast] and associated  diffusion
%                  coefficients for ground truth (combined in fdInput[])
% or
% - dRange[] and = reasonable intervals for the random choice of f and d
%   fRange[]       made by DiffParamRandomizer.m
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Init parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
InitVar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inputSimu = zeros(6,length(b),numberOfIter);
fitNNLS = zeros(4,m,numberOfIter);

for i = 1:numberOfIter
    % Randomize diffusion parameters gT in specified intervals every iter
    [d,f, fdInput, DBasis, rawSignal] = createSignal(d,f,b , Dmin, Dmax, m);

    % Running NNLS simulations (w & w\ reg) generating new noise every iter
    [signal, noise, snr]       = createNoise(rawSignal, snrIn); 

    [sNNLSNoReg, sNNLSReg, mu] = NNLSfitting(DBasis, signal, i/numberOfIter*100);
    inputSimu(:,:,i) = [b; fdInput; noise(:)'; signal'; snr']; 
    fitNNLS(:,:,i)   = [DValues; sNNLSNoReg'; sNNLSReg'; [mu zeros(1,m-1)]];
end

% Calculating NNLS diffusion parmeters (0 = noReg, 1 = Reg)
[dNNLS, fNNLS, resultNNLS, dAUC, fAUC] = findpeaksNNLS(fitNNLS, np, 1);

% NLLS* with NNLS results as a priori information
 [dNLLS, fNLLS, resnormNLLS]  = NLLSfitting(inputSimu, Dmin, Dmax, dNNLS, fNNLS);
 resultNLLS(:,:)              = [dNLLS; fNLLS; resnormNLLS];

% NLLS with standard starting values
[dNLLSraw, fNLLSraw, resnormNLLSraw]  = NLLSfitting(inputSimu, Dmin, Dmax);
resultNLLSraw(:,:)                    = [dNLLSraw; fNLLSraw; resnormNLLSraw];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting and saving data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting figures
plotSimu(inputSimu, fitNNLS, dNNLS, dNLLS, fNLLS, numberOfIter);
[simResults,statistic] = evaluation(resultNNLS, resultNLLS, resultNLLSraw, dAUC, fAUC, fdInput, np, snrIn, m, fd_offset);

% Write simulation data to file
write3DMatrixToTxt(inputSimu, "results/SimuInput.txt" ); 
write3DMatrixToTxt(fitNNLS,   "results/NNLSfit.txt"   ); 
writetable        (simResults,"results/simResults.xlsx"); 