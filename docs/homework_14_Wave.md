# **计算物理第14次作业　Wave**  
`作者：陈洋遥`  `学号2013301020169`  `更新时间：201605031` 

## **内容目录**  
[TOC]
## **本文摘要** 
　　波动现象是自然界十分常见的现象，从弹性波(例如绳波、水波等)到电磁波，再到引力波，波动现象是物理学研究的重要课题。解决波动问题，首先要建立波动方程，而后的问题就是在给定边界条件下求解波动方程。本次作业完成课后习题6.6,6.13，展示绳波的传播、反射过程，给出给定点振动的功率谱。
　
## **研究背景** 
　　绳波是一维波动问题，在不考虑能量损耗、认为绳是轻柔绳的前提下，绳波由如下波动方程描述：
\begin{equation}
\frac{\partial^2 y}{\partial t^2}=c^2 \frac{\partial^2 y}{\partial x^2}
\tag{1}
\end{equation}
其中$c$是波速。要解出绳波的运动状态，必须给定一定的边界条件。常见的边界条件有自由边界条件和固定边界条件等。为简单起见，这里取固定边界条件，即绳端总是不发生横向位移。对于非端点的元段，则可以让其离散化，然后通过如下的迭代方法逐步求解绳波随时间的演化：
　　　　$y(i,n+1)=2[1-r^2]y(i,n)$
\begin{equation}
-y(i,n-1)+r^2[y(i+1,n)+y(i-1,n)]
\tag{2}
\end{equation}
这里$r=c\Delta t/\Delta x$，$t=n\Delta t$和$x=i \Delta x$分别是时间和元段位置坐标。初始条件通过如下方法给出：在连续两个时刻，给出绳波位移，这样即可开始迭代，给出绳波随时间的演化。
## **绳波问题的求解** 
### **绳波随时间的演化**  
　　以一个高斯波包来描述绳波被激发时的形态，通过求解上面的波动方程，给出绳波的演化如下：
<center>![中心激发高斯波包](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter6_wave/ch3_wave_middle.gif) </center><center>`图1 中心激发高斯波包`</center><center>　　图片展示在绳子中心激发一个高斯波包后绳波随时间的变化。可以发波包一分为二，反向运动，在固定边界发生反射后往回传播，再次相遇。</center> 

<center>![两个高斯波包 ](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter6_wave/ch6_wave_twopeak.gif) </center><center>`图2 两个高斯波包`</center><center>　　图片展示同时激发两个高斯波包的情况。可以发现，两个波包均分裂成子波包，在各波包相遇时，振动叠加，而在分离后，各自以原来的形态传播，波包互相不发生影响。这是由波动方程的线性性决定的。</center> 

<center>![正弦波激发](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter6_wave/ch6_wave_sine.gif) </center><center>`图3 正弦波激发`</center><center>　　图片展示以正弦波的形式激发绳波。可以发现，绳波直接是一个驻波。</center> 

<center>![正弦波包](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter6_wave/ch6_wave_packet.gif) </center><center>`图4 正弦+高斯波包激发`</center><center>　　如图展示一个高斯波包组合上正弦波形成的一个复合波包的演化。可以发现，波包仍然是逐渐一分为二，遇到固定边界时发生反射。</center> 

<center>![将绳上一点拉开](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter6_wave/ch3_wave_triangle.gif) </center><center>`图5 将绳上一点拉起`</center><center>　　图片展示将绳上一点拉起后，绳的形态随时间的变化。可以发现，绳子在演化过程中能够保持折线状态。</center> 
### **绳上给定点位移的功率谱** 
　　将绳上某点的振动信号作快速傅里叶分解并取模平方，得到信号的功率谱展示如下。在不同位置激发、在不同位置观察时，功率谱呈现不同特征，有的峰位可以缺失。
<center>![ ](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter6_wave/ch6_signal_03_01.png) </center><center>`图6 在离端点10%位置观察绳波`</center><center>　　图片展示在离端点10%位置观察绳波的功率谱，绳子用高斯波包激发，位置在离端点30%位置，可以发现5倍频和10倍频都是缺失的</center> 
<center>![ ](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter6_wave/ch6_signal_03_02.png) </center><center>`在离端点20%位置观察绳波`</center><center>　　　　图片展示在离端点20%位置观察绳波的功率谱，绳子用高斯波包激发，位置在离端点30%位置，可以发现，仍然是5倍频和10倍频缺失</center> 
<center>![ ](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter6_wave/ch6_signal_05_05.png) </center><center>`在离端点50%位置观察绳波`</center><center>　　　　图片展示在离端点50%位置观察绳波的功率谱，绳子用高斯波包激发，位置在离端点50%位置。可以发现，只能激发奇数倍频，且5倍频还是缺失。</center> 
<center>![ ](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter6_wave/ch6_signal_gauss_sine.png) </center><center>`图9 高斯+正弦波包激发的频谱`</center><center>　　图片展示高斯+正弦波包激发的绳波在给定观察点的功率谱</center> 

## **小结**
　　本次作业展示绳波的运动过程。观察了不同激发形式的绳波随时间的演化过程。给出了在绳上给定点观察绳波信号的功率谱。

## **致谢和引用**
[1] 计算物理；Nicholas J. Giordano, Hisao Nakanishi.
[2] 常用数学符号的LaTex表示方法；[http://www.mohu.org/info/symbols/symbols.htm](http://www.mohu.org/info/symbols/symbols.htm).
[3] matplotlib-绘制精美的图表；[http://old.sebug.net/paper/books/scipydoc/matplotlib_intro.html](http://old.sebug.net/paper/books/scipydoc/matplotlib_intro.html).


























