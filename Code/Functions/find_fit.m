function [y_corrected,ym,z,F,gof,fitted_curve]=find_fit(f,channel,check,ch,ML)
% Use this function to find the fit for the 1/f component of the power spectrum. The output variables are:
% y_corrected: the power spectrum after removing the 1/f component
epochs_p=channel(ch).power;
F=repmat(f',[size(epochs_p,1),1]);
fi=(f>3).*(f<250);
fi=logical(fi);
if ML~=1
    check=~check;
end
Y=epochs_p(check,fi);
X=F(check,fi);
x=double(X(:));
y=double(Y(:));
% x0=[1,1e5,2];
fitfun = fittype( @(b,k,a,x) b-log10(k+x.^a));
yy=log10(1e6*y);
[fitted_curve,gof] = fit(x,yy,fitfun,'Robust','BiSquare','Algorithm','Trust-Region',Lower=[-1e3,-1e3,0]);
Z=fitted_curve(X(1,:));
z=10.^(Z)/1e6;
ym=mean(Y,1);
y_corrected=mean(Y,1)-z';
F=X(1,:);
coeffvals = coeffvalues(fitted_curve);
format bank
disp(gof.rsquare)
disp(coeffvals)