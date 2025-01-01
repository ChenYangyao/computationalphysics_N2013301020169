# **计算物理第11次作业　The Solar System**  
`作者：陈洋遥`  `学号2013301020169`  `更新时间：20160508` 


## **内容目录**  
[TOC]
## **本文摘要**  
　　本次作业研究太阳系中行星运动问题。相比于地球上的物体运动受到显著的空气阻力影响，行星运动是良好的物理观测场所。17世纪开普勒发现行星运动三大定律，而后牛顿又成功地用万有引力定律解释了开普勒定律，行星运动的研究推动了物理学的发展。本次作业完成课后习题4.7,4.9；顺带给出广义相对论引力理论的一个实验验证实例。
## **研究背景**  
　　牛顿的万有引力定律给出了两个质点之间的引力相互作用的规律。在太阳系中，行星受到太阳的引力作用为   
\begin{equation}
f_G=\frac{GM_S M_P}{r^2}
\tag{1}
\end{equation}此式连同牛顿第二定律，就可以给出行星运动规律。倘若考虑太阳的有限质量，则行星与太阳二体问题也容易解决，只需要把惯性质量$M_E$全部换成所谓的约化质量$\mu=M_SM_P/(M_S+M_P)$，则仍可形式上地利用牛顿第二定律解决行星相对太阳的运动问题。  
　　通过简单的积分可以得到行星运动的轨道方程为  
\begin{equation}
r=\frac{l}{1-e\cos\theta}
\tag{2}
\end{equation}其中$l=L^2/(\mu G M_P M_S)$,$e$是轨道离心率，其取值可以决定轨道性态。由这样的轨道方程就可以给出开普勒的三大定律。需要注意的是，开普勒三定律是牛顿引力理论即平方反比引力的结果，若力不是平方反比的，则开普勒定律就不成立，行星也一般不存在闭合轨道。下面就来逐步研究行星运动问题。
## **双星问题**  
　　前面已经提到双星问题可以合理的化为单星问题，但星体数量更多时就不能化为简单的运动问题。因此有必要在静止坐标系下直接对双星运动求数值解。我们展示求解结果如下图1所示，相应的程序见[ch4_solar_sys_binary.py](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter4_solar_sys/ch4_solar_sys_binary.py)。可以发现，双星运动轨迹高度依赖初始条件。
<center>![双星问题](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter4_solar_sys/ch3_binary.png) </center><center>`图1 牛顿引力下的双星问题`</center><center>图展示牛顿引力理论下的双星问题，图中两星体质量满足$M_1/M_2=2$，以不同的初始条件开始运动。左图(a)显示初始条件合适时，双星轨道为圆形，右图(b)显示初始条件偏离成圆条件时，双星轨道为椭圆。计算中单位均取天文学单位。</center>    

　　我们也展示双星运动的动态过程如下，其参数取值和图1中的(b)图是类似的，双星均沿椭圆轨道运动<center>![双星动态图](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter4_solar_sys/ch3_binary.gif) </center><center>`图2 双星问题动态图 `</center><center></center>  

## **引力不是平方反比的情况**  
　　如果引力不是按牛顿引力理论预言的那样是平方反比的，而是另外的比例系数，即$f_G=GM_SM_P/r^\beta$，$\beta\neq2$，这样来讲一般行星轨道将产生进动。我们利用数值方法在不同初始条件下求解此时的行星运动，相应的程序见[ch4_solar_sys_inversesquare.py](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter4_solar_sys/ch4_solar_sys_inversesquare.py)。计算结果如下图3所示。
<center>![非平方反比引力下星体的运动](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter4_solar_sys/ch3_non_square_precession.png) </center><center>`图3 非平方反比引力下星体的运动 `</center><center>图展示引力系数不是平方反比情况的星体运动，这里引力是$\beta=2.05$次方反比的，认为中心天体不动。图(a)展示了初始条件合适时，行星仍然能够保持圆轨道运动，此时不发生进动。而图(b)(c)则展示了行星轨道偏离圆周运动时，轨道将不是闭合的椭圆，而是将不断进动，离心率越大时，进动越明显。图(d)给出了进动速率随离心率变化的规律。</center>  

　　我们也展示上面过程的动态图形如下，其中参数取值同上面图3(a)(b)(c)三图。<center>![非平方反比引力下星体的运动](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter4_solar_sys/precession_pass.gif) </center><center>`图4 非平方反比引力下星体的运动动态图 `</center><center>图中很明显地展示了行星轨道在非平方反比引力作用下的进动，且展示了离心率与进动的关系。</center>  
　　
## **相对论引力效应**
　　众所周知，行星运动的牛顿理论在很大程度上是成功的，几乎可以解释当时所观测到的大部分天文现象。然而水星的进动的理论却一直难以和实验观测符合，即便考虑所有能对水星运动产生显著影响的因素，理论和实验仍然有每世纪43分的偏差。然而，爱因斯坦的引力理论却成功解释了水星的进动。然而这里我们要展示相对论引力的另外的两个效应——光线的偏折(1919，Eddington)和引力导致的红移(1959，Pound & Rebka)，这两个现象完全无法用牛顿引力理论解释。  
　　根据爱因斯坦广义相对论引力理论，球对称星体周围的稳态真空度规由Schwarzchild解给出，在自然单位制下，度规形式为  
\begin{equation}
\mathrm{d}s^2=-(1-\frac{2M}{r})\mathrm{d}t^2+(1-\frac{2M}{r})^{-1}\mathrm{d}r^2+r^2\mathrm{d}\Omega^2
\tag{2}
\end{equation}Schwarzchild度规完全确定了星体周围的时空结构，而光子在这样的时空中将沿着类光测地线运动，其运动方程由测地线方程给出   
\begin{equation}
\frac{\mathrm{d}^2x^\mu}{\mathrm{d}\lambda^2}+\Gamma^{\mu}_{\sigma\rho}\frac{\mathrm{d}x^\sigma}{\mathrm{d}\lambda}\frac{\mathrm{d}x^\rho}{\mathrm{d}\lambda}=0  
\tag{2}
\end{equation}其中$\Gamma^{\mu}_{\sigma\rho}$是Christoffel Symbols，上式连同类光测地线切失的归一化条件即可确定光线运动的路径。而离星球不同距离的稳态观察者发射和接收光信号的频率是不一样的，若接收者离星体较远，则其接收到的光信号频率将变低，具体来讲为$\omega_{source}/\omega_{receiver}=(1-2M/r_{receiver})^{1/2}/(1-2M/r_{source})^{1/2}$。下面就利用数值方法求解光信号的运动路径，并大致地展示光型号的红移过程。结果如下图5所示。
<center>![光信号的偏折和红移](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter4_solar_sys/redshift_pass.gif) </center><center>`图5 大质量恒星附近光信号的偏折和红移 `</center><center>如图所示，光源(Source)靠近大质量星体，其发出的光信号将在弯曲时空中产生偏折，而离星体较远的接受者(Receiver)接收到的信号频率将降低，即发生了红移。光的偏折和红移纯粹是相对论引力效应。</center> 


    
## **小结**
　　本次作业讨论了太阳系中行星运动问题。通过求解两体运动方程得到了两体运动在静止坐标系下的解，通过研究非平方反比引力，展示了行星进动一个可能的来源。最后展示了非线性引力对光信号的偏折和红移效应。

## **致谢和引用**
[1] 计算物理；Nicholas J. Giordano, Hisao Nakanishi.
[2] 常用数学符号的LaTex表示方法；[http://www.mohu.org/info/symbols/symbols.htm](http://www.mohu.org/info/symbols/symbols.htm).
[3] matplotlib-绘制精美的图表；[http://old.sebug.net/paper/books/scipydoc/matplotlib_intro.html](http://old.sebug.net/paper/books/scipydoc/matplotlib_intro.html).












