function obj = horzcat(adObj, varargin)
% HORZCAT(varargin) concatenates the property arrays of each ADvar object
% in VARARGIN to the current object ADOBJ.

% TODO: why bother defining this as a class function if it does not modify
% the ADvar object itself?

    obj = ADvar(adObj.val, adObj.der);
    
    for k = 1:length(varargin)
        list_element = varargin{k};
    
        if ~isa(list_element, 'ADvar')
            obj.val = horzcat(obj.val, list_element);
            obj.der = horzcat(obj.der, zeros(size(list_element)));
    
        else
            obj.val = horzcat(obj.val, list_element.val);
            obj.der = horzcat(obj.der, list_element.der);
    
        end
    end
    
end