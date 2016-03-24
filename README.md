# filtfilt
We wirte a filtflit function in java . The filtflit's output is the same as it's in Matlab .


最近实验室做android项目需要使用零相位滤波，在Matlab上有直接的实现，
C++也有相应的版本，唯独java中在网上没有找到。因而萌生的自己动手写
一个java实现的想法。

本项目参考了在网上搜的一个C++对于filtfilt的实现，在文件夹“filtfilt-C++”中有visual studio 2010的工程
由于原作者不可考证，在这里特别感谢这位匿名的同学。


在文件夹“iltfilt-Java”中有自己对filtfilt的java实现。

-----------------------------测试代码-----------------------------

		ArrayList<Double> B = new ArrayList<Double>();
		ArrayList<Double> A = new ArrayList<Double>();
		ArrayList<Double> X = new ArrayList<Double>();
		
		B.add(4.0);
		B.add(5.0);
		B.add(6.0);
		
		A.add(1.0);
		A.add(2.0);
		A.add(3.0);
		
		X.add(1.0);
		X.add(2.0);
		X.add(3.0);
		X.add(4.0);
		X.add(5.0);
		X.add(6.0);
		X.add(7.0);
		X.add(8.0);
		
		ArrayList<Double> y = Filtfilt.doFiltfilt(B,A,X);
		for (int i = 0; i < y.size(); i++)
			System.out.println(y.get(i));
			
	---------------------------输出结果---------------------------
				-6731884.25
				7501778.75
				-2757230.25
				-662443.25
				1360955.75
				-686678.25
				4135.75
				227147.75                              与Matlab的输出完全相同。




