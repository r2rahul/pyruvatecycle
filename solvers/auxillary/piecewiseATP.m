function out = piecewiseATP(glc, atp)
%Function for modeling ATP dynamics
if glc < 10e-3
out = atp;
elseif glc >= 10e-3
out = 7e-3;
else
out = 7e-3;
end
