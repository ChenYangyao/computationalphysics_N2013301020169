# **计算物理第7次作业　Baseball: Motion of a Batted Ball**  
`作者：陈洋遥`  `学号2013301020169`  `更新时间：20160408`  
## **内容目录**  
[TOC]
## **本文摘要**  
　　本次作业编写程序解决棒球运动问题。棒球运动的特点是：球速高，球轨迹变化多。棒球在空中飞行时，其空气拖拽系数将严重依赖于速度，并且棒球旋转将使其受到额外的Magnus力，因而棒球的轨迹可以变得十分复杂。  
　　本次作业主要解决课本2.19习题，附带讨论：  
　　(1)空气阻力对棒球轨迹、射程的影响；
　　(2)弧线球的成因和其对棒球旋转角速度的依赖。

## **研究背景**  
　　棒球的飞行是典型的抛体运动，但是棒球在飞行过程中除受到重力作用外，还会受到空气阻力和因棒球旋转产生的Magnus力，一般来讲，棒球飞行轨迹由下述方程给出
\begin{equation}
m \frac{\mathrm{d} ^2 \vec r}{\mathrm{d}t ^2}=m \vec g - B_2v^2 \vec v^0+S_0 \vec v \times \vec \omega,
\tag{1}
\end{equation}  
其中$m$是棒球质量，$g$是重力加速度，$\vec v^0$是速度方向的单位矢量，$\omega$即棒球旋转角速度。上面$-B_2 v^2 \vec v^0$即空气阻力项，而$S_0 \vec v \times \vec \omega$即为Magnus力。在棒球速度范围内尚可认为$S_0$基本是一常量，而$B_2$则显著依赖于棒球的速度，按照教科书上所给近似，可将其写为  
\begin{equation}
\frac{B_2}{m}=0.0039+\frac{0.0058}{1+e^{\frac{v-35}{5}}}.(\mathrm{SI})
\tag{2}
\end{equation}  
式$(1)(2)$结合初始条件即可解得棒球飞行轨迹，但却难以给出解析解。不妨利用欧拉法求其近似解，可给出其相应的差分方程为  
\begin{equation}
\left\{
\begin{aligned}  
x_{i+1}&=v_{x,i}\Delta t+x_{i}  \\  
v_{x,i+1}&=(-(B_2/m) v v_x+(S_0/m) v_y \omega _z) \Delta t+v_{x,i}  \\
y_{i+1}&=v_{y,i}\Delta t+y_{i}  \\
v_{y,i+1}&=(-g-(B_2/m) v v_y+(S_0/m) v_z \omega _x) \Delta t+v_{y,i} \\
z_{i+1}&=v_{z,i} \Delta t+z_i \\
v_{z,i+1}&=(S_0/m) v_x \omega _y \Delta t+v_{z,i}  \\
\end{aligned}
\right. 
\tag{3}
\end{equation}  
　　只要知道棒球前一时刻的状态，就可通过上述方程给出下一时刻的状态，从而反复迭代最终获得棒球的整个运动轨迹。由于上面的方程是近似成立的，因而不可避免的要引入误差。欧拉法每一步的误差量级是$\Delta t^2$，而在给定总时间内的总体误差将是$\Delta t$量级，因此适当缩小步长$\Delta t$将有助于提高运算精度，但往往会使计算时间以$\frac{1}{\Delta t}$增加，因此需要权衡两者利弊，取一在给定时间不产生显著误差的的步长即可，下面的计算将会体现这一点。  
## **棒球运动的研究**  
* **空气阻力对棒球运动的影响**
　　为研究空气阻力对棒球运动的影响，我们只需分别计算有、无阻力时棒球运动的情况即可。一般来说无阻力的棒球应该能飞得更远，下面的计算验证了这一点。在方程$(3)$中我们取棒球初速度为$110 \mathrm{mph}$，发射角$\theta = 45 ^o$ ,阻力系数$B_2/m$由式$(2)$给出，而Magnus力的系数则取为$S_0/m=4.1 \times 10^{-4}$,对无空气阻力的情况只需把两个系数均取为0即可，程序见[ch2_baseball_backspin.py](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter2_201604_12/ch2_baseball_backspin.py)。我们展示计算结果如下图由图所示。为确定步长$\mathrm{d} t$的合理取值，我们分别用不同的步长计算了抛物结果，当步长逐渐减小时，结果逐渐收敛到一稳定结果，可以发现步长取$0.10\mathrm{s}$和$0.05 \mathrm{s}$时棒球轨迹不再明显变化，因此后面的计算中均取$\mathrm{d} t= 0.10 \mathrm{s}$。从图中可以发现:  
  - 空气阻力大大减少了棒球的射程，有空气阻力时棒球轨迹明显低于无空气阻力的情况，另外，棒球的旋转对于其射程也有不小的影响。  
  - 空气阻力的存在也很大程度上影响棒球的射程，通过具体计算(见程序[ch2_baseball_range_comp.py](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter2_201604_12/ch2_baseball_range_comp.py))可知其射程在有、无阻力时分别见下表。可见空气阻力将棒球射程减小了一半还要多。
  |有无阻力| 发射角|最大射程|
|:----:| :---:| :--:|
|无空气阻力|$45 ^o$|$252 \mathrm{m}$|
|有空气阻力|$35 ^o$|$123 \mathrm{m}$|  
![棒球的运动](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter2_201604_12/ch2_fig_batted_ball.png)  
　　　　　　　　　　　　　　　　　　`图1  有、无空气阻力时棒球的运动 `   

*  **棒球在3d空间的运动**  　　
　　下面研究棒球的3d运动，基本参数取值与前面相同，只不过这次棒球旋转的角速度不再沿水平方向，而是沿竖直方向，因此棒球的轨迹将显著偏离发射面，此时棒球将在3d空间运动。同理可用欧拉法解决方程$(3)$，相应的程序见[ch2_baseball_backspin_3d.py](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter2_201604_12/ch2_baseball_backspin_3d.py)，棒球运动偏离发射面的程度由其旋转的角速度决定，可见角速度达到$400 \mathrm{rad/s}$时偏转甚至可以达到$50 \mathrm{m}$以上，棒球的旋转显著影响其运动轨迹。
　　　　　　![棒球在3d空间的运动](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter2_201604_12/ch2_fig_batted_ball_3d.png)  
　　　　　　　　　　　　　　　　　　　　`图2  棒球在3d空间的运动 `　　
　　　　　　　　`棒球旋转角速度沿竖直方向，因而运动偏离发射面，且偏离由角速度大小决定`　　
　　
## **小结**
　　本次作业讨论了棒球的运动问题，讨论了棒球所受空气阻力和Magnus力对棒球运动的影响。其中前者的影响极为显著，将很大地减少棒球的射程，明显地改变棒球的轨迹，而后者也有不小的影响，如果棒球旋转适当，可以朝不同的方向发生偏移，因而给棒球运动带来很大的不确定性和娱乐性。

## **致谢和引用**
[1] 计算物理；Nicholas J. Giordano, Hisao Nakanishi
[2] 常用数学符号的LaTex表示方法；[http://www.mohu.org/info/symbols/symbols.htm](http://www.mohu.org/info/symbols/symbols.htm)
[3] matplotlib-绘制精美的图表；[http://old.sebug.net/paper/books/scipydoc/matplotlib_intro.html](http://old.sebug.net/paper/books/scipydoc/matplotlib_intro.html)
[4] 参考了刘星辰大神的3d作图方法，这里表示$10000^{100}$分感谢，顺祝好人一生平安




