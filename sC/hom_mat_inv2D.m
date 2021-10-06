function GstInv = hom_mat_inv2D(Gst)
           %
           %
           %
            Rst = Gst(1:2,1:2);
            dst = Gst(1:2,3);
            GstInv = [Rst.',-Rst.'*dst;0,0,1];
            
            
        end