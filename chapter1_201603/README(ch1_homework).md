# 计算物理Chapter1(第4次) 作业  
## 摘要  
　　本次作业内容为教材课后习题1.5：两种原子核各有几率衰变成另一种核，因此它们的衰变规律将由两个互相耦合的一阶线性微分方程描述。我采用课本上介绍的欧拉方法来解决此问题，并用图像表示出了两种核数目随时间的变化规律。  
　　另外，根据直观感觉，欧拉法的计算精度应该随着步长的减小而增加，因此计算中采用了几种不同的步长，并通过图形展示出了不同步长下求得的数值解与解析解的偏差，计算发现，当步长小于衰变常数的1/20时，欧拉法的计算精度是十分可观的。  
　　通过本次作业的练习，我终于基本掌握了用python绘图的方法，另外，也第一次通过自己写程序来解决微分方程。    
## 背景  
　　课本习题1.5 Consider again  a decay problem with two types of nuclei A and B, but now suppose that nuclei of type A decay into ones of type B, while nuclei of type B decay into ones of type A. Strictly speaking, this is not a "decay" process, since it is possible for the type B nuclei to turn back into type A nuclei. A better analogy would be a resonance in which a system can tunnel or move back and forth between two states A and B which have equal energies. The correponding rate equation are   
　　　　　　　　　　　　　　　d N_A /dt = N_B/tau - N_A/tau  　　
　　　　　　　　　　　　　　　　　d N_B /dt = N_A/tau - N_B/tau　    
　
where for simplicity we have assumed that the two types of decay are characterized by the same time constant, tau. Solve this system of equations for the numbers of nuclei, N_A and N_B, as functions of time. Consider different initial conditions, such as N_A=100, N_B=0, etc., and take tau =1s. Show that you numerical results are consistent with the idea that the system reaches a steady state in which N_A and N_B are constant. In such a steady state, the time dericatives dN_A/dt and dN_B/dt should vanish.  
 

## 正文  
　　源代码（[戳我戳我](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/homework/chapter1_201603/chapter1_homework_20160327_1.py)），计算得到的数据([pickle形式](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter1_201603/nuclei_decay_two_type_txt.txt)，[txt数据](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter1_201603/nuclei_decay_two_type_txt.txt))，截图([指令窗](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter1_201603/ch1_2.png),[数据结果](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/homework/chapter1_201603/decay%20of%20nuclei.jpg))    
 　　为解决课本上的习题，需用到的的函数库有  
　　　　　　`import pickle`；  
　　　　　　`from pylab import *`；  
　　　　　　`import numpy as py`；  
　　其中第一个是用来导入、导出数据的；第二个则是绘图的函数库；第三个用来制造数组极为方便。  
　　欧拉法的原理是将微分方程化为差分方程，并利用前一时刻的结果来计算后一时刻的结果，反复迭代即得到变量对时间的依赖关系。欧拉法要能求得比较精确的解，步长的选取是十分关键的。一般来说，选取原则是：步长必须小于系统中任意的特征参数(例如，本题中就是要小于衰变的特征时间)。  
　　结果分析：本次计算中，衰变特征时间tau=10，取步长分别为四种情况：5.0， 3.0， 2.0， 0.5，  由此分别计算出数值解并作图(如下左图)。另外，为更直观地表示出各情况数值计算的精度，分别计算了各情况数值解与解析结果的偏差(如右下图)。可以发现，在步长为时间常数的1/2时,数值解偏离正确值很大，在核数为1000左右的情况下，计算偏差能达到200；而当步长减小时，计算误差则也减小；当步长仅为时间常数的1/20时，数值偏差衰减到不足20，此时可以认为计算精度比较满意。    
　!(衰变图片)[https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/chapter1_201603/decay%20of%20nuclei.jpg]  
## 致谢  
　　感谢caihao老师的作业例子(点击进入)[https://github.com/caihao/computational_physics_whu/tree/master/chapter1]  
　　matplotlib的教程(点击进入)[http://liam0205.me/2014/09/11/matplotlib-tutorial-zh-cn/]
　　python的数据读取和写入教程(点击进入)[https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/homework/chapter1_201603/decay%20of%20nuclei.jpg]
=======
# 计算物理Chapter1(第4次) 作业  
## 摘要  
　　本次作业内容为教材课后习题1.5：两种原子核各有几率衰变成另一种核，因此它们的衰变规律将由两个互相耦合的一阶线性微分方程描述。我采用课本上介绍的欧拉方法来解决此问题，并用图像表示出了两种核数目随时间的变化规律。  
　　另外，根据直观感觉，欧拉法的计算精度应该随着步长的减小而增加，因此计算中采用了几种不同的步长，并通过图形展示出了不同步长下求得的数值解与解析解的偏差，计算发现，当步长小于衰变常数的1/20时，欧拉法的计算精度是十分可观的。  
　　通过本次作业的练习，我终于基本掌握了用python绘图的方法，另外，也第一次通过自己写程序来解决微分方程。  
## 背景  
　　课本习题1.5 Consider again  a decay problem with two types of nuclei A and B, but now suppose that nuclei of type A decay into ones of type B, while nuclei of type B decay into ones of type A. Strictly speaking, this is not a "decay" process, since it is possible for the type B nuclei to turn back into type A nuclei. A better analogy would be a resonance in which a system can tunnel or move back and forth between two states A and B which have equal energies. The correponding rate equation are   
　　　　　　　　　　　　　　　d N_A /dt = N_B/tau - N_A/tau  
　　　　　　　　　　　　　　　d N_B /dt = N_A/tau - N_B/tau  
where for simplicity we have assumed that the two types of decay are characterized by the same time constant, tau. Solve this system of equations for the numbers of nuclei, N_A and N_B, as functions of time. Consider different initial conditions, such as N_A=100, N_B=0, etc., and take tau =1s. Show that you numerical results are consistent with the idea that the system reaches a steady state in which N_A and N_B are constant. In such a steady state, the time dericatives dN_A/dt and dN_B/dt should vanish.
## 正文  
　　源代码（[戳我戳我](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter1_201603/chapter1_homework_20160316_2.py)），计算得到的数据([pickle形式](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter1_201603/nuclei_decay_two_type_pickle.txt)，[txt数据](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter1_201603/nuclei_decay_two_type_txt.txt))，截图([指令窗](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter1_201603/ch1_2.png),[数据结果](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/chapter1_201603/decay%20of%20nuclei.jpg))  
　　在第二次修改中，给图增加了单位，调整了字体使之更规范。  
 　　为解决课本上的习题，需用到的的函数库有  
　　　　　　`import pickle`；  
　　　　　　`from pylab import *`；  
　　　　　　`import numpy as py`；  
　　其中第一个是用来导入、导出数据的；第二个则是绘图的函数库；第三个用来制造数组极为方便。  
　　欧拉法的原理是将微分方程化为差分方程，并利用前一时刻的结果来计算后一时刻的结果，反复迭代即得到变量对时间的依赖关系。欧拉法要能求得比较精确的解，步长的选取是十分关键的。一般来说，选取原则是：步长必须小于系统中任意的特征参数(例如，本题中就是要小于衰变的特征时间)。  
　　本次计算中，衰变特征时间tau=10，取步长分别为四种情况：5.0, 2.0， 1.0， 2.0，0.5, 由此分别计算出数值解并作图(如下左图)。另外，为更直观地表示出各情况数值计算的精度，分别计算了各情况数值解与解析结果的偏差(如右下图)。可以发现，在步长为时间常数的1/2时,数值解偏离正确值很大，在核数为1000左右的情况下，计算偏差能达到200；而当步长减小时，计算误差则也减小；当步长仅为时间常数的1/20时，数值偏差衰减到不足20，此时可以认为计算精度比较满意。  
　　![衰变图片](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/homework/chapter1_201603/decay%20of%20nuclei.jpg) 
## 致谢  
　　感谢caihao老师的作业例子[点击进入](https://github.com/caihao/computational_physics_whu/tree/master/chapter1)  
　　matplotlib的教程[点击进入](http://liam0205.me/2014/09/11/matplotlib-tutorial-zh-cn/)  
　　python的数据读取和写入教程[点击进入](http://www.ibm.com/developerworks/cn/opensource/os-python8/)
>>>>>>> origin/master
