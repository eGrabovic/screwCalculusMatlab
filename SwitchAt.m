function res = SwitchAt(v1, v2, k)
% SWITCHAT(v1, v2, k) is an internal utility function.
%   TODO: what does it do?
		res = v1;
		index = k;
		res(index) = v2(index);
end