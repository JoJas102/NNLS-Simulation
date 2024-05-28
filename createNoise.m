function [signal, noise, SNRout] = createNoise(signalIn, snrIn)

    % Generating noise (as overlay) new for every iteration
    noiseRaw = randn(size(signalIn));                 % noise based on gaussian distribution of randn()
    noise_factor = signalIn(1)/(snrIn*std(noiseRaw));
    noise = noise_factor*noiseRaw;
    signal = signalIn + noise_factor*noiseRaw;
    SNRout = signal./abs(std(noise_factor*noiseRaw)); % consider first b-value [Periquito2021]
    %signal = signalIn; fprintf('ATTENTION: no noise is added to signal \n'); % un-comment for run w/o noise
end