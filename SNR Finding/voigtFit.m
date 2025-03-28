% Fit a voigt function to data using a minimum log likelihood fit. "data"
% is the data that needs fitting, xx is the corresponding x vector.

function [peak,fwhm,chi2] = voigtFit(data,xx)
    
    initial = [2e9 1e9 xx(round(length(xx)/2)) 0.13];
    param = fminsearch(@(var) sum((data-voigt(var(1),var(2),var(3),var(4),xx)).^2),initial);
    fitVoigt = voigt(param(1),param(2),param(3),param(4),xx);
    chi2 = sum((data-fitVoigt).^2./fitVoigt);
    
    numPts = round((xx(end)-xx(1))/1e6);
    x = linspace(xx(1),xx(end),numPts);
    fitVoigt = voigt(param(1),param(2),param(3),param(4),x);
    peak = max(fitVoigt);
    
    temp = fitVoigt > peak/2;
    strt = find(temp,1);
    final = find(temp,1,'last');
    fwhm = x(final)-x(strt);
    
end

