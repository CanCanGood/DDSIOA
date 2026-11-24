# Case 1：运行与绘图

使用 Windows 和 MATLAB R2024a，将 MATLAB 当前文件夹切换到本目录。

## 运行

```matlab
result = run_case1('GFDHS', 0.05, 1, 'results/gfdhs_run1');
```

四个输入依次为：算法名称、以米表示的空间设置、随机种子、结果保存目录。结果目录必须是尚不存在的新目录。

算法可选 `GFDHS`、`GBH`、`Astar`。空间设置可选 `0.01`、`0.05`、`0.10`、`0.15`、`0.20`、`0.25`。

```matlab
run_case1('GBH',   0.25, 1, 'results/gbh_run1');
run_case1('Astar', 0.25, 1, 'results/astar_run1');
```

入口会自动读取 `data/` 中的数据并选择对应算法函数，不要将所有子目录同时加入 MATLAB 路径。建议先运行以上示例；细网格运行需要更多内存。

## 绘图

```matlab
plot_case1('results/gfdhs_run1/result.mat');
```

运行结果保存在指定目录中，包括 MAT/JSON 文件和路径坐标 CSV 文件，导出坐标的单位为米。
