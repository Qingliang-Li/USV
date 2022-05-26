# USV
此程序对应文章：李清亮, 李彬, 孙国皓, 等. 基于精确罚函数的无人艇航迹规划和自动避障算法 [J]. 中国舰船研究, 2021, 16(1): 89–95.
LI Q L, LI B, SUN G H, et al. Trajectory planning and automatic obstacle avoidance algorithm for unmanned surface vehicle based on exact penalty function[J]. Chinese Journal of Ship Research, 2021, 16(1): 89–95.
说明：该问题采用精确罚函数处理约束，使用控制参数化方法求解最优控制问题

miser3是程序启动的主函数；
fig是画图的函数；
ocf:写动态方程
ocg0:写终端约束
ocg:写不等式约束
ocphi:终端等式约束
ocxzero:初始值

使用说明：
运行miser3.m;
勾选“Read in and execute an existing data input file”;
键入文本名，如“T1”;
等待计算结果
