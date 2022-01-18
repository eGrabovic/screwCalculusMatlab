function res = TensorDerivative(J, q)
    DJ  = diff(J, q);
    res =(1/2) * (DJ + DJ');
end