% fast vectorized implementation for sampling from shifted wald distributions
% all parameters are vectors of equal size. For details refer to:
% Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating random variates
% using transformations with multiple roots. The American Statistician, 30(2), 88-90.
% (c)-2012 Amiya Patanaik
function samples = shiftedWald(Ter,mu,lambda)
    %sample from a normal distribution with a mean of 0 and 1 standard deviation
    sampleSize = length(mu); 
    v = randn(sampleSize,1);
    y = v.*v;
    samples = mu + (mu.*mu.*y)./(2*lambda) - (mu./(2*lambda)) .* sqrt(4*mu.*lambda.*y + mu.*mu.*y.*y);
    test = rand(sampleSize,1);
    alternate = (mu.*mu)./samples;
    altIndex = test > mu./(mu+samples);
    samples(altIndex) = alternate(altIndex);
    samples = samples + Ter;
end