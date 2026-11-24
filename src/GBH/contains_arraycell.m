function result = contains_arraycell(y, target)
target=target(1:3) ; 
len = length(target);
  [rows, cols] = size(y);
  result = false;
  
  for i = 1:rows
    for j = 1:cols-len+1                              
      if y(i, j:j+len-1) == target
        result = true;
        return;
      end
    end
  end
end