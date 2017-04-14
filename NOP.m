% This function returns an array that is laid out as follows: row number is
% which element we are in (following the order of iterate in z, then y,
% then x). Columns 1 through 8 are the global node numbers of the nodes
% that make up that element
function [NOP_array] = NOP(NNY, NNZ, NNX)
num_elements = (NNY-1)*(NNZ-1)*(NNX-1);
NOP_array = zeros(num_elements, 8);
element =1;
for elx=1:NNX-1
    for ely=1:NNY-1
        for elz=1:NNZ-1
            NOP_array(element, 1) = (elx-1)*NNY*NNZ + (elz-1)*NNY + ely;
            NOP_array(element,2) = (elx-1)*NNY*NNZ + (elz-1)*NNY + ely + 1;
            NOP_array(element,3) = (elx-1)*NNY*NNZ + (elz-1)*NNY +ely + NNY;
            NOP_array(element,4) = (elx-1)*NNY*NNZ + (elz-1)*NNY +ely + NNY + 1;
            NOP_array(element,5) = (elx-1)*NNY*NNZ + (elz-1)*NNY + NNZ*NNY + ely;
            NOP_array(element,6) = NOP_array(element,5) + 1;
            NOP_array(element,7) = NOP_array(element,5) + NNY;
            NOP_array(element,8) = NOP_array(element,6) + NNY;
            
            element= element+1;
        end
    end
end

end