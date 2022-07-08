function [AEP_MWh, Energy] = AEP(WSP,P,M,k)
    A                  = M/(gamma(1+1/k));
    CDF                = 1-exp(-(WSP/A).^k);
    N_CDF              = length(CDF);
    for ii = 1:N_CDF-1
        Time(ii)        = CDF(ii+1)-CDF(ii);
        MeanP(ii)       = (P(ii+1)+P(ii))/2;
    end
    Energy             = Time.*MeanP;
    AEP_MWh            = sum(Energy)*8760/1000;
end