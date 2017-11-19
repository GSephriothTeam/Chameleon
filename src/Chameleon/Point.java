package Chameleon;

import java.util.ArrayList;
import java.util.List;

//package DataMining_Chameleon;

public class Point{
	//坐标点id号,id号唯一
	int id;
	//point value
	List<Integer> list; 
	
	//坐标横坐标
	//int x;
	//坐标纵坐标
	//int y;
	//是否已经被访问过
	boolean isVisited;
	
	public Point(String s){
		//this.id = Integer.parseInt(id);
		//this.x = Integer.parseInt(x);
		//this.y = Integer.parseInt(y);
		String []ss = s.split(" ");
		this.id = Integer.parseInt(ss[0]);
		List<Integer> l = new ArrayList<Integer>();
		for ( int i = 1; i < ss.length; i ++)
			l.add(Integer.parseInt(ss[i]));
		this.list = l;
	}
	
	/**
	 * 计算当前点与制定点之间的欧式距离
	 * 
	 * @param p
	 *            待计算聚类的p点
	 * @return
	 */
	public double ouDistance(Point p) {
		double distance = 0;
		
		for ( int i = 0; i < p.list.size() ; i ++)
			distance += (this.list.get(i) - p.list.get(i)) * (this.list.get(i) - p.list.get(i));
		//distance = (this.x - p.x) * (this.x - p.x) + (this.y - p.y) * (this.y - p.y);
		distance = Math.sqrt(distance);

		return distance;
	}
	
	/**
	 * 判断2个坐标点是否为用个坐标点
	 * 
	 * @param p
	 *            待比较坐标点
	 * @return
	 */
	public boolean isTheSame(Point p) {
		return this.list.equals(p.list);
		/*
		boolean isSamed = false;

		if (this.list.equals(p.list)) {
			isSamed = true;
		}
		return isSamed;
		*/
	}
}
