function v = TakeMoments(A)
% Parameters reshaping
		v = [A(1,1), -A(1,2), -A(1,3), A(2,2), -A(2,3), A(3,3)];
end
