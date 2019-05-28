function Aav = calc_frechet_mean( A, wgt, verbose )
%CALC_FRECHET_MEAN  find the Frechét mean of a set of transformation matrices
% 
%   Aav = CALC_FRECHET_MEAN( A ) returns average transformation A of a
%   cell vector of transformations A
%
%   CALC_FRECHET_MEAN( A, wgt ) applies weights in numeric vector wgt
% 
%   CALC_FRECHET_MEAN( A, wgt, verbose ) outputs verbose messages if
%   verbose is true

% jfpva (joshua.vanamerom@kcl.ac.uk)


%% Info

% Compute the average transformation
% Please refer to the following papers for more information
% [1] Marc Alexa. Linear combination of transformations.
% [2] Maria Kuklisova-Murgasova et al. A dynamic 4D probabilistic atlas of the developing brain.
% [3] Paul Aljabar et al. Assessment of brain growth in early childhood using deformation-based morphometry.


%% Initialise

if ~exist( 'wgt', 'var' )
    wgt = ones( size( A ) );
end

if ~exist( 'verbose', 'var' )
    verbose = false;
end

normLogDeltaMu = inf;
tolerance = 0.0000000000001;
iterations = 20;

totalWeight = sum( wgt );

if ( totalWeight <= 0 )
    error('Sum of weight must be positive.');
end

mu = A{1};

n = 0;


%% Find Mean

while ( normLogDeltaMu > tolerance && n < iterations )

    sumLogs = zeros(4);

    for i = 1:numel(A)
        deltaM = mu \ A{i};
        sumLogs = sumLogs + wgt(i) * logm( deltaM );
    end

    sumLogs = sumLogs / totalWeight;
     
    deltaMu = expm( sumLogs );

    mu = mu * deltaMu;

    logDeltaMu = logm( deltaMu );
    
    normLogDeltaMu = sum(abs(logDeltaMu(:)));
    
    if verbose
        fprintf( 'IRTK FrechetMean (iteration %i, normLogDeltaMu=%g)\n', n, normLogDeltaMu )
        disp(mu)
    end

    n = n + 1;

end

Aav = mu;

if verbose
    fprintf( 'IRTK FrechetMean\n' )
    disp(Aav)
end


end  % calc_frechet_mean(...)