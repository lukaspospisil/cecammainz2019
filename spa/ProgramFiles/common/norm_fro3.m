function [ aa ] = norm_fro3( A )
%NORM_FRO3 Summary of this function goes here
%   Detailed explanation goes here

aa = 0;
for i=1:size(A,3)
    aa = aa + norm(A(:,:,i),'fro')^2;
end

aa = sqrt(aa);

end

