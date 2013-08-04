function [P, data] = diis(data, P, S, X, F)
	if isempty(data)
        data = struct('P', {}, 'e', {});
    end
	data(end+1).P = P;
	data(end).e = X'*(F*P*S-S*P*F)*X;
	if length(data) > 5
		data(1) = [];
	end
	while true
		B = form_B({data.e});
		s = svd(B);
		if s(end)/s(1) > 1e-15 % check if invertible
            break
        end
		data(1) = [];
	end
	m = length(data);
	c = B\[zeros(m, 1); -1]; % obtain coeffs
	P = zeros(size(data(1).P));
	for i = 1:m
		P = P + c(i)*data(i).P; % DIIS trial density
	end
end

% forms B matrix from [Pulay1980]
function B = form_B(e)
	m = length(e);
	B = zeros(m);
	for i = 1:m
		for j = 1:i
			bij = sum(sum(e{i}.*e{j}));
			B(i,j) = bij;
			B(j,i) = bij;
		end
	end
	B = [B -ones(m, 1); -ones(1, m) 0];
end
