function h_val = h(cur_node, goal)
h_val = min(sum(abs(repmat(cur_node,size(goal,1),1)-goal),2 ));

end