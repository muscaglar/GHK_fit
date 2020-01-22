    function [output] = size_check(in1,in2,in3,in4)
       if(size(in1,1)==1)
       else
           in1 = in1';
       end
       
       if(size(in2,1)==1)
       else
           in2 = in2';
       end
       
       if(size(in3,1)==1)
       else
           in3 = in3';
       end
       
       if(size(in4,1)==1)
       else
           in4 = in4';
       end
       output = [in1',in2',in3',in4'];
    end