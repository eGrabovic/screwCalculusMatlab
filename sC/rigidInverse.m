function GstInv = rigidInverse(Gst)
            %
            % Computes the inverse of a homogeneous transormation Gst into
            % Gst^-1
            %
%             if all(size(Gst) ~= [4,4])
%                 error('Homogenous transformation matrix requires to have 4x4 dimension')
%             end
if isa(Gst, 'ADvar')
    
    Der = rigidInverse(Gst.der);
    GstInv = ADvar(rigidInverse(Gst.val), Der);
    return
    
end
            Rst = Gst(1:3,1:3);
            dst = Gst(1:3,4);
            GstInv = [Rst.',-Rst.'*dst; 0 0 0 Gst(4,4)];
end