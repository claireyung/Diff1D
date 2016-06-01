%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     DIFF1Dredout              %
% This script reduces the size of the output of %
% Diff1Dsolve.m                                 %
%                                               %
% Ryan Holmes March 2016                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (NOUT ~= 1)
%Create averaging matrix:
NtOUT = floor(Nt/NOUT);
avMAT = zeros(NtOUT,Nt);
for ii = 1:NtOUT
    avMAT(ii,(1+(ii-1)*NOUT):(ii*NOUT)) = 1/NOUT;
end
avMAT = avMAT';

%reduce size of variables:
vars = {'t','shflux','srflux','ssflux','tau_x',...
    'tau_y','Hsbl','u','v','w','S','T','b','ks','kt','kv',...
    'Tadv','gams','gamt','gamv','dbdy','bulkRiN','bulkRiD'};
for vi = 1:length(vars)
    if (exist(vars{vi}))
        eval([vars{vi} ' = ' vars{vi} '*avMAT;']);
    end
end

clear avMAT vars;

end
