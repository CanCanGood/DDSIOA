function random_array = generate_sum_to_one_array(n)
  % 生成n-1个随机数
  rng('default');
  random_numbers = rand(1, n);
  % 归一化随机数，使其和为1
  random_array = random_numbers / sum(random_numbers);

end