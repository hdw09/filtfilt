# filtfilt
We wirte a filtflit function in java . The filtflit's output is the same as it's in Matlab .


���ʵ������android��Ŀ��Ҫʹ������λ�˲�����Matlab����ֱ�ӵ�ʵ�֣�
C++Ҳ����Ӧ�İ汾��Ψ��java��������û���ҵ�������������Լ�����д
һ��javaʵ�ֵ��뷨��

����Ŀ�ο����������ѵ�һ��C++����filtfilt��ʵ�֣����ļ��С�filtfilt-C++������visual studio 2010�Ĺ���
����ԭ���߲��ɿ�֤���������ر��л��λ������ͬѧ��


���ļ��С�iltfilt-Java�������Լ���filtfilt��javaʵ�֡�

-----------------------------���Դ���-----------------------------

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
			
	---------------------------������---------------------------
				-6731884.25
				7501778.75
				-2757230.25
				-662443.25
				1360955.75
				-686678.25
				4135.75
				227147.75                              ��Matlab�������ȫ��ͬ��




