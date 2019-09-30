function [ normA ] = gershgorin( A )

normA_row = diag(A) - abs(diag(A)) + sum(abs(A),2);
normA = max(normA_row);

end

