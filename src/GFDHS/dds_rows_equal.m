function equal=dds_rows_equal(rows,point)
% Exact single-row membership, preserving NaN inequality and input order.
if isempty(rows),equal=false(size(rows,1),1);return;end
equal=all(rows==point,2);
end
