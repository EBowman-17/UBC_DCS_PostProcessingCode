function voigtFunc = voigt(sigma,gamma,peakCenter,A,xx)
    %VOIGT, this function creates a voigt profile based on the sigma and gamma
    %passed in as input arguments. xx is an x vector for the data.

    %First calculate a z vector for the voigt function
    z = ((xx-peakCenter)+gamma*1i)/(sqrt(2)*sigma);
    voigt = real(fadf(z))/sigma;
    voigtFunc = A*voigt/max(voigt);

end

