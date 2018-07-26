function SNR_mat = calcSNR(F_mat,timeDim,MODEFLAG,baseline)
%This function calculates a SNR matrix with time along columns. User must
%input the time dimension of input matrix (timeDim) and the baseline period
%over which to average F0 unless MODEFLAG is set to 1.
%
%SNR matrix is calculated as DFF/sqrt(F0), which should serve as a valid
%measure except for such cases when there is a dark current value that is
%close to F0 (i.e., N~F0 so that DFF*N^(1/2)~DFF*F0^(1/2)).
%
%If MODEFLAG is set to 1, values within the bottom 25 percentile are
%averaged for F0. This is useful in cases where heightened fluorescence is
%sparse as is typical for, e.g. spontaneous astrocyte calcium data.

if timeDim~=1 %if time dose not lie along columns, take matrix transpose
    F_mat = F_mat';
end

if ~MODEFLAG
    %Calc DFF based on declared baseline SD
    F0_vect = mean(F_mat(1:baseline,:),1);
    F0_mat = repmat(F0_vect,[size(F_mat,1),1]);
    DFF_mat = (F_mat - F0_mat)./F0_mat;
    DFF_std = std(DFF_mat(1:baseline,:),1);
    DFF_stdMat = repmat(DFF_std,[size(F_mat,1),1]);
    SNR_mat = DFF_mat./DFF_stdMat;
    
elseif MODEFLAG
    %Get upper and lower limits to average for DFF
    SNR_mat = zeros(size(F_mat));
    for iROI = 1:size(F_mat,2)
        F = F_mat(:,iROI);
        F_bQuart = quantile(F,0.25);
        F0 = mean(F(F<F_bQuart));
        if F0==0
            F0 = 0.1;
        end
        SNR_mat(:,iROI) = (F - F0)/sqrt(F0);
    end
    
end