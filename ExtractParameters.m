function p = ExtractParameters(CGtable, Masslist, Tensortable) 
	p = [];
	n = size(CGtable,1);
	
    for k = [1, n]
		
         p0k = [Masslist(k)];
         
         p1k = Masslist(k) * CGtable(k);
         
         p2k = TakeMoments(Tensortable(k) -...
             Masslist(k) * Hat(CGtable(k)) * Hat(CGtable(k)));
         				
		 pk = [p0k, p1k, p2k];
		 
		 p = [p, pk]; 
    end 
end