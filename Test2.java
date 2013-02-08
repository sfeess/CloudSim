import java.awt.BorderLayout;
import java.awt.Graphics;

import javax.swing.*;



public class Test2  {
	
	

	 public static void main(String[] args) {
		 float test;
		 
		 test= -5.123456789123456789f;
		 test -= (int) (test);
		 test=test*(float) Math.pow(10, 7);
		 test = (int)test;
		 System.out.println(test);
	   }
	 }