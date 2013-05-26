import java.awt.Color;



public class Test2  {
	
	

	 public static void main(String[] args) {
		 
		 	int ssx=10;
		 	int ssy=10;
		 
		 	float xpos=0.9f;
		 	float ypos=5;
		 	
		 	
		 	int x1 = (int) xpos;
		 	int x2 = x1+1;
		 	int y1 = (int) ypos;
		 	int y2 = y1+1;
			
		 	System.out.println("x1="+x1+"   x2="+x2+"   y1="+y1+"   y2="+y2);
		 	
			// Periodic X-Boarder
			if(x1<1){					x1 = ssx + x1;	
				if(x2<1)				x2 = ssx + x2;	
			}
			
			else if(x2>ssx){			x2 = (x2-ssx);	
				if(x1>ssx)				x1 = (x1-ssx);	
			}			
			
			
			System.out.println("x1="+x1+"   x2="+x2+"   y1="+y1+"   y2="+y2);
			
			// Y-Boarder
			if(y1<0){					y1 = 0;	}
			if(y2<0){					y2 = 0;}
			if(y1>ssy+1){				y1 = ssy+1;}			
			if(y2>ssy+1){				y2 = ssy+1;}
			
			System.out.println("x1="+x1+"   x2="+x2+"   y1="+y1+"   y2="+y2);
			
			// Distances
			float dx = (float) (xpos-Math.floor(xpos));//xpos%1;
			float dy = (float) (ypos-Math.floor(ypos));//ypos%1;
		 
		 
			
		 
		 
		 
		 
		 
	   }
	 }