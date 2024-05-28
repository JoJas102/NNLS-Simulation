    format long

    % Initialize variables
    numberOfIter = 1000;   % number of simulations
    snrIn = 140;
    b = [0,5,10,20,30,40,50,75,100,150,200,250,300,400,525,750]; % [Stabinska2020] + low b-values vital [LeBihan2017]

    n = length(b); 
    f_ref = [60 30 10];        % 60 30 10 % fill up with zeroes when simulating less components
    d_ref = [1 5.8 165].*1e-3; % 1 5.8 165
    f = f_ref;
    d = d_ref;
    dRange = [0.7 1.5;     % reasonable intervals for the choice of d and f
                3 7;       
              70 200].*1e-3;
    fRange = [55 75;       
              15 45;
               3 15];    
    np = nnz(d);           %np = length(dRange); % number of input peaks

    % Randomize diffusion parameters in specified intervals
    %[d,f] = gTRandomiser(dRange, fRange); fprintf('ATTENTION: you are fitting wiht random ground truth values! \n');

    % Generate basis values
    Dmin = 0.7*1e-3;       % standard [vanBaalen2016] & [Wong2019]
    Dmax = 300*1e-3;
    m = 350;               % no. of bins             
    DValues = logspace(log10(Dmin), log10(Dmax), m); 
    DBasis = exp(-kron( b', DValues));
    
    % Real d-values (exact bins)
    [dMin(1), index(1)] = min(abs(DValues-d(1)));
    [dMin(2), index(2)] = min(abs(DValues-d(2)));
    [dMin(3), index(3)] = min(abs(DValues-d(3)));
    d = [DValues(index(1)) DValues(index(2)) DValues(index(3))];
    fdInput = [f zeros(1,n-length(f)); d zeros(1,n-length(d))];
    
    % Calculate offset to reference values for evaluation
    fd_offset = [d - d_ref 0;
                 f - f_ref 0];

    % Define total diffusion signal decay according to sum of f_i*e^(-b*D_i)
    A = exp(-kron(b', d)); % constraint matrix containing exp decay funxtions
    rawSignal = A*f';      % mean signal of ROI/voxel without noise


    % Variations

    %b = [0,5,10,20,30,40,50,75,100,150,200,250,300,400,525,750]; % [Stabinska2020] + low b-values vital [LeBihan2017]
    %b = [0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750]; % equidistant
    %b = [0,5,10,15,20,30,50,100,150,200,250,350,450,550,650,750]; % block/interval
    %b = [0,5,10,15,20,25,35,50,70,100,140,200,275,384,535,750]; % fast focus (=interpretation of exp. sequence exp(linspace(log(5),log(750),16))
    %b = [0,5,15,30,40,50,60,70,80,90,100,150,200,300,500,750]; % inter focus
    %b = [0,50,105,160,210,265,320,370,425,480,430,575,635,690,745,750]; % slow focus (=reverse exp)
    %b = [0,5,10,20,30,50,75,100,150,200,250,300,450,600,750,900]; % medium extended b-range [Periquito2020]
    %b = [0,5,10,20,30,50,75,100,150,200,300,400,550,700,850,1000]; % extended b-range [Periquito2020]
    %b = [0,5,10,20,30,50,75,100,150,250,350,500,650,800,1000,1200]; % extreme extended b-range [Periquito2020]

    %Dmin = 0.9*1e-3; % shortened extreme
    %Dmax = 200*1e-3;
    %Dmin = 0.5*1e-3; % extended extreme
    %Dmax = 500*1e-3;