package com.iiphci.dsp;

import java.util.ArrayList;

/**
 * 零相位滤波函数
 * 
 * @author 黄大卫,胡盼
 * @Email hhhdddwww@163.com
 * @version 1.0
 * @Time 2016年3月23日下午9:04:28
 * @Description 调用static方法：doFiltfilt执行滤波
 * @copyright 安徽大学 iiphci 实验室
 * 
 */

public class Filtfilt {

	/**
	 * 
	 * @param B
	 *            滤波器参数B
	 * @param A
	 *            滤波器参数A
	 * @param input
	 *            要滤波的数组 ： 要求是一维，长度要大于B数组。
	 * @return 滤波后的结果 返回double数组。
	 */
	public static ArrayList<Double> doFiltfilt(ArrayList<Double> B,
			ArrayList<Double> A, ArrayList<Double> X) {

		int len = X.size();
		int na = A.size();
		int nb = B.size();
		int nfilt = (nb > na) ? nb : na;
		int nfact = 3 * (nfilt - 1);
		if (len <= nfact)
			throw new RuntimeException("输入数值X长度太小，数据最少是滤波器阶数的三倍");
		resize(B,nfilt, 0);
		resize(A,nfilt, 0);
		
		ArrayList<Integer> rows = new ArrayList<Integer>();
		ArrayList<Integer> cols = new ArrayList<Integer>();

		add_index_range(rows, 0, nfilt - 2,1);
		if (nfilt > 2)
		{
			add_index_range(rows, 1, nfilt - 2, 1);
			add_index_range(rows, 0, nfilt - 3, 1);
		}
		add_index_const(cols, 0, nfilt - 1);
		if (nfilt > 2)
		{
			add_index_range(cols, 1, nfilt - 2,1);
			add_index_range(cols, 1, nfilt - 2,1);
		}
		int klen = rows.size();
		ArrayList<Double> data = new ArrayList<Double>();
		resize(data,klen, 0);
		data.set(0, 1 + A.get(1));
		int j = 1;
		if (nfilt > 2)
		{
			for (int i = 2; i < nfilt; i++)
				data.set(j++, A.get(i));
			for (int i = 0; i < nfilt - 2; i++)
				data.set(j++, 1.0);
			for (int i = 0; i < nfilt - 2; i++)
				data.set(j++, -1.0);
		}
		ArrayList<Double> leftpad = subvector_reverse(X, nfact, 1);
		changeArray2(leftpad,2*X.get(0));
		
		ArrayList<Double> rightpad = subvector_reverse(X, len - 2, len - nfact - 1);
		changeArray2(rightpad,2*X.get(len - 1));
		
		double y0;
		ArrayList<Double> signal1 = new ArrayList<Double>();
		ArrayList<Double> signal2 = new ArrayList<Double>();
		ArrayList<Double> zi = new ArrayList<Double>();
		//reserve(signal1,leftpad.size() + X.size() + rightpad.size(),0);
		append_vector(signal1, leftpad);
		append_vector(signal1, X);
		append_vector(signal1, rightpad);
		
		double [][] sp = Zeros(max_val(rows)+1,max_val(cols)+1); 
		for (int k = 0; k < klen; ++k)
			sp[rows.get(k)][cols.get(k)] = data.get(k);
		double[]bb = map(B);
		double[]aa = map(A);
		//Mat.inv(sp)
		//bb.segment(1, nfilt - 1) - (bb(0) * aa.segment(1, nfilt - 1))
		
		double[][] ZZi = Mat.multi(Mat.inv(sp),calc(segment(bb,1,nfilt - 1),bb[0],segment(aa,1,nfilt - 1)));
		
		resize(zi,ZZi.length,1);
		
		changeZi(ZZi,zi,signal1.get(0));
		filter(B, A, signal1, signal2, zi);
		reverse(signal2);
		changeZi(ZZi,zi,signal2.get(0));
		filter(B, A, signal2, signal1, zi);
		ArrayList<Double> Y = subvector_reverse(signal1, signal1.size() - nfact - 1, nfact);
		return Y;

	}

	private static void reverse(ArrayList<Double> signal2) {
		int i=0;
		int j=signal2.size()-1;
		while(i<j){
			swap(signal2,i,j);
			i++;
			j--;
		}
	}

	private static void swap(ArrayList<Double> signal2, int i, int j) {
		double temp = signal2.get(j);
		signal2.set(j, signal2.get(i));
		signal2.set(i, temp);
	}

	private static void changeZi(double[][] zZi, ArrayList<Double> zi,
			Double double1) {
		for (int i = 0; i < zZi.length; i++) {
			zi.set(i, zZi[i][0]*double1);
		}
		
	}

	private static double[][] calc(double[] segment, double d, double[] segment2) {
		double[][] ret = new double[segment.length][1];
		for (int i = 0; i < segment.length; i++) {
			ret[i][0] = segment[i]-d*segment2[i];
		}
		return ret;
	}

	private static double[] segment(double[] bb, int i, int j) {
		double[]ret=new double[j-i+1];
		for (int k = 0; k < j-i+1; k++) {
			ret[k] = bb[i+k];
		}
		return ret;
	}

	private static double[] map(ArrayList<Double> b) {
		double[] ret = new double[b.size()];
		for (int i = 0; i < ret.length; i++) {
			ret[i] = b.get(i);
		}
		return ret;
	}

	private static double[][] Zeros(int ii, int jj) {
		double [][] sp = new double[ii][jj];
		for (int i = 0; i < ii; i++)
			for (int j = 0; j < jj; j++)
				sp[i][j] = 0;
		return sp;
	}

	public static void filter(ArrayList<Double> B, ArrayList<Double> A,
			ArrayList<Double> X, ArrayList<Double> Y,ArrayList<Double> Zi) {
		if (A.size() == 0)
			throw new RuntimeException("A 数组为空！");
		boolean flagA = true;
		for (Double doubleA : A) {
			if (doubleA != 0) {
				flagA = false;
			}
		}
		if (flagA) {
			throw new RuntimeException("A 数组至少要有一个数不为零！");
		}
		if (A.get(0) == 0) {
			throw new RuntimeException("A 数组第一个元素不能为零！");
		}
		changeArray(A, A.get(0));
		changeArray(B, A.get(0));
		
		
		int input_size = X.size();
		int filter_order = max(A.size(),B.size());
		resize(B,filter_order,0);
		resize(A,filter_order,0);
		resize(Zi,filter_order,0);
		resize(Y,input_size,0);
		
		for (int i = 0; i < input_size; i++) {
			int order = filter_order - 1;
			while(order!=0){
				if(i>=order)
					Zi.set(order-1, B.get(order)*X.get(i-order)-A.get(order)*Y.get(i-order)+Zi.get(order));
				--order;
			}
			Y.set(i, B.get(0)*X.get(i)+Zi.get(0));
		}
		Zi.remove(Zi.size()-1);		
	}

	

	private static void resize(ArrayList<Double> a, int i, double j) {
		if(a.size()>=i)
			return;
		int size = a.size();
		for (int j2 = size; j2 < i; j2++) {
			a.add(j);
		}
	}

	private static int max(int size, int size2) {
		if(size >size2)
			return size;
		else
			return size2;
	}

	static void changeArray(ArrayList<Double> vec, double a0) {
		for (int i = 0; i < vec.size(); i++) {
			vec.set(i, vec.get(i)/a0);
		}
	}
	
	static void changeArray2(ArrayList<Double> vec, double a0) {
		for (int i = 0; i < vec.size(); i++) {
			vec.set(i, a0-vec.get(i));
		}
	}

	static void add_index_range(ArrayList<Integer> indices, int beg, int end,
			int inc) {
		for (int i = beg; i <= end; i += inc)
			indices.add(i);
	}

	static void add_index_const(ArrayList<Integer> indices, int value, int numel) {
		while (numel-- != 0)
			indices.add(value);
	}

	static void append_vector(ArrayList<Double> vec, ArrayList<Double> tail) {
		for (Double doubleitem : tail) {
			vec.add(doubleitem);
		}
	}

	static ArrayList<Double> subvector_reverse(ArrayList<Double> vec,
			int idx_end, int idx_start) {
		ArrayList<Double> resultArrayList = new ArrayList<Double>(idx_end
				- idx_start + 1);
		for (int i = 0; i < idx_end- idx_start + 1; i++) {
			resultArrayList.add(0.0);
		}
		int endindex = idx_end - idx_start;
		for (int i = idx_start; i <= idx_end; i++)
			resultArrayList.set(endindex--, vec.get(i));
		return resultArrayList;
	}

	static int max_val(ArrayList<Integer> vec) {
		int temp = vec.get(0);
		for (Integer integer : vec) {
			if (temp < integer)
				temp = integer;
		}
		return temp;
	}

}
