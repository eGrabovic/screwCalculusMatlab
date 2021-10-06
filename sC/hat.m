        function X = hat(vec)
            %
            % trasforma un vettore colonna R3 in forma matriciale hat R3x3
            % antisimmetrica
            %
            % oppure trasforma un twist colonna R3 in forma omogenea hat
            % R4x4
            %         

                if isa(vec, 'ADvar')
                    X = ADvar(hat(vec.val), hat(vec.der));
                    return
                end
                
                if size(vec, 1) == 3
                    X = [0 -vec(3) vec(2);...
                        vec(3) 0 -vec(1);...
                        -vec(2) vec(1) 0];
                    return
                end
                
                X = [hat([vec(4);vec(5);vec(6)]), vec(1:3);0 0 0 0];


        end