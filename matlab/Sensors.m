function [BB, pqr] = Sensors(BB, pqr)
    for i = 1:3
        %%% MagField Sesor
        MagBias = 4e-9; % Tesla
        MagFieldBias = MagBias*(2*rand() - 1);
        MagNoise = 1e-7; % Tesla
        MagFieldNoise = MagNoise*(2*rand() - 1);

        %%% Angle Rates Sesor
        AngleScaleBias = 9*pi/180/3600; % 9deg/hour
        AngleBias = AngleScaleBias*(2*rand() - 1);
        AngleScaleNoise = 0.15*pi/180/sqrt(3600); % 0.15deg/sqrt(hour)
        AngleNoise = AngleScaleNoise*(2*rand() - 1);
        
        BB(i) = BB(i) + MagFieldBias + MagFieldNoise;
        pqr(i) = pqr(i) + AngleBias + AngleNoise;
    end
    
end