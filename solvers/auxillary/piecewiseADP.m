function out = piecewiseADP(glc, adp)
%Function for modeling ADP dynamics
if glc < 10e-3
out = adp;
elseif glc >= 10e-3
out = 0.6e-3;
else
out = 0.6e-3;
end