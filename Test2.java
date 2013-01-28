import java.awt.BorderLayout;
import java.awt.Graphics;

import javax.swing.*;



public class Test2  {
	
	
	
	static int ssx = 10000;
	static int ssy = 10000;
	static float x;
	
	static float[] onedim = new float[ssx*ssy];
	static float[][] twodim = new float[ssx][ssy];
	
	public static int flin(int i, int j){		
		return ((i)+(ssx)*(j));}
	
	 public static void main(String[] args) {
		long start, time;
		System.out.println("Comparison of 1D and 2D Arrays (size=10000x10000)");
		//********************************************************
		start = System.nanoTime();
		 for(int i= 0; i<ssx; i++){
			for(int j= 0; j<ssy; j++){
				onedim[flin(i,j)] = 176.5349f+i+j;
			}
		}
		time = System.nanoTime() - start;
		System.out.println("write 1D: "  + time/1000000f+ " ms");

		//********************************************************
		start = System.nanoTime();
		 for(int i= 0; i<ssx; i++){
				for(int j= 0; j<ssy; j++){
					twodim[i][j] = 176.5349f+i+j;
				}
			}
		time = System.nanoTime() - start;
		System.out.println("write 2D: "  + time/1000000f+ " ms");
		
		//********************************************************
		start = System.nanoTime();
		for(int i= 0; i<ssx; i++){
				for(int j= 0; j<ssy; j++){
					x= onedim[flin(i,j)];
					//System.out.println(onedim[flin(i,j)]);
				}
			}
		time = System.nanoTime() - start;
		System.out.println("read 1D: "  + time/1000000f+ " ms");
 
		//******************************************************** 
		start = System.nanoTime();
		for(int i= 0; i<ssx; i++){
				for(int j= 0; j<ssy; j++){
					x= twodim[i][j];
					//System.out.println(twodim[i][j]);
				}
			}
		time = System.nanoTime() - start;
		System.out.println("read 2D: "  + time/1000000f+ " ms");
 
		 
		 
	   }
	 }