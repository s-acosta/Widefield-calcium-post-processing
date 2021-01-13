function getMemSize(this)
    
   props = properties(this); 
   totSize = 0; 
   
   for ii=1:length(props) 
      currentProperty = getfield(this, char(props(ii))); 
      s = whos('currentProperty'); 
      totSize = totSize + s.bytes; 
   end
  
   totSize = totSize / 1e9;
   fprintf(1, '%.2f GB\n', totSize); 

end