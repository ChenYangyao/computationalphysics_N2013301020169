# 计算物理 第三次作业
>2013301020169 陈洋遥
20160309 更新  

## 摘要  
本文展示了第三次作业内容，实现了以下功能：  

* 在屏幕上用字符串组成的阵列输出自己姓名的功能，且输出顺序可以是任意的；  
* 在屏幕上以字符串阵列的形式画出了一条正弦曲线，且可以随时间推移向左移动

## 背景  


## 正文  
* 在屏幕上用字符串阵列的形式显示自己的姓名。实现此功能基本原理是利用Python内的字典命令`{key:value}`把字母用字符串阵列的形式存于计算机内，然后按字母检索相应阵列，并用`print`指令一行一行输出该阵列
    * 代码在此([戳我戳我](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/homework20160309_1.py))  
    * 结果如图，用户输入字符chen yy后屏幕上即显示出相应的字符串阵列，而若打乱输入顺序，则仍然能够成功输出；
    ![字符串输出图](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/03.png)
    * 值的改进的地方：若输入字符太多，超过一行的显示范围，则不能成功显示，可以增加自动换行功能；另外，该程序只能显示少数几个字母以及空格，若要显示其他的字母，则需要一一添加至字典；  

* 在屏幕上用字符阵列画出自己想画出的东西，并可以让"图形"动起来。其基本原理和显示字母是类似的，将图形以"像素点"的形式储存于字符阵列即可。另外，根据老师的提示，使用指令`import os`，`i = os.system('cls')`即可清空屏幕上已经显示的"图形"，再输出另外的"帧",即可让图形动起来。再者，使用指令`import time`，`time.sleep(sometime)`即可让计算机输出每帧的间隔拉长，避免动画太快。
    * 代码在此([戳我戳我](https://github.com/ChenYangyao/computationalphysics_N2013301020169/blob/master/homework20160309_2.py))  
    * 结果如图，成功在屏幕上输出一条"正弦曲线"(至少看起来有点像)，并且，该曲线可以随着时间推移逐渐向左移动(这里只把其中一"帧"截图放了上来)。该图形模拟了一束向-x方向传输的正弦波。
    ![正弦波传输图](https://raw.githubusercontent.com/ChenYangyao/computationalphysics_N2013301020169/master/03.2.png)  
    * 值的改进的地方：该图形的点过于分散，不够密，因为每个字符本来就有那么大，所以暂时没想到改善的方法。  

## 结论  
本次作业实验了利用字符串阵列输出字母或简单图形的功能，但编程过程较复杂，期待利用python的画图功能进一步改善。

## 致谢  
* 感谢caihao老师关于清空屏幕的指令
`import os`，`i = os.system('cls')`
的作业提示([戳我戳我](https://github.com/caihao/computational_physics_whu/blob/master/Exercises.md))  
* 感谢百度知道提供让python暂停工作的指令([点击戳戳百度知道](http://zhidao.baidu.com/link?url=CKjyQhhvlTis0udu69dltkFlvIwtHtJRPE91z4gQzPl6jhbCuE03aZQ9gVGQB7P4eBiyAtapdsvg7HDOpH0GEq))
`import time`，`time.sleep(sometime)`




