
public class Interpolation {

	static int ssx,ssy;
	static int x1,x2,y1,y2;			// interpolation positions
	static float dx,dy;				// interpolation distances
	static float xpos, ypos;		// interpolation results

	
	public static void setSpace(int x, int y){
		ssx=x;
		ssy=y;
	}
	
	public static float biLin(float xpos, float ypos, float[][] f){
		
		x1 = (int) xpos;
		x2 = x1+1;
		y1 = (int) ypos;
		y2 = y1+1;
		
		if(x2>ssx+1) {x1 = ssx+1; 	x2 = ssx+1;	}
		if(y2>ssy+1) {y1 = ssy+1; 	y2 = ssy+1;	}
		if(x1<0)	 {x1 = 0; 		x2 = 0;		}
		if(y1<0)	 {y1 = 0; 		y2 = 0;		}
		
		dx = xpos%1;
		dy = ypos%1;
		
		return  	(1-dx) *( 	f[x1][y1]*(1-dy) 	+ 	f[x1][y2]*dy  	)
					+  dx  *(	f[x2][y1]*(1-dy)	+  	f[x2][y2]*dy  	);
	}
	
	
	private static float cubic(float x, float[] f){
		float d_i, d_i1, d_q;
		
		d_i	 = (f[2] - f[0])	/2;
		d_i1 = (f[3] - f[1])  	/2;
		d_q  = (f[2] - f[1]);
		
		//if signs different -> slopes to zero
		if ((d_q >= 0 && d_i < 0)  || ((d_q< 0 && d_i >= 0)))  	d_i	 =0;
		if ((d_q >= 0 && d_i1 < 0) || ((d_q< 0 && d_i1 >= 0)))  d_i1 =0;
		
		return  f[1] + d_i*x + (3*d_q - 2*d_i-d_i1)*x*x + (-2*d_q + d_i + d_i1)*x*x*x;
	}
	
	
	public static float biCubic(float xpos, float ypos, float[][] f){
		int x = (int)xpos;
		int y = (int)ypos;
		
		int x1 = x-1;
		int x2 = x;
		int x3 = x+1;
		int x4 = x+2;
		
		int y1 = y-1;
		int y2 = y;
		int y3 = y+1;
		int y4 = y+2;
		
		
		if(x<1){x1=0;}
		if(y<1){y1=0;}
		
		if(x<0){x2=0; x3=1; x4=2;}
		if(y<0){y2=0; y3=1; y4=2;}
		
		if(x>ssx-1){x4=ssx+1;}
		if(y>ssy-1){y4=ssy+1;}
		
		if(x>ssx){x3=ssx+1;}
		if(y>ssy){y3=ssy+1;}
		
		if(x>ssx+1){x2=ssx+1; x1=ssx;}
		if(y>ssy+1){y2=ssy+1; y1=ssy;}
		
		
		float[] arr_x = new float[4];
		float[] arr_y = new float[4];
		
		//interpolate first line
		arr_x[0] = f[x1][y1];
		arr_x[1] = f[x2][y1];
		arr_x[2] = f[x3][y1];
		arr_x[3] = f[x4][y1];
		arr_y[0] = cubic(xpos%1, arr_x);
		
		//interpolate second line
		arr_x[0] = f[x1][y2];
		arr_x[1] = f[x2][y2];
		arr_x[2] = f[x3][y2];
		arr_x[3] = f[x4][y2];
		arr_y[1] = cubic(xpos%1, arr_x);
		
		//interpolate third line
		arr_x[0] = f[x1][y3];
		arr_x[1] = f[x2][y3];
		arr_x[2] = f[x3][y3];
		arr_x[3] = f[x4][y3];
		arr_y[2] = cubic(xpos%1, arr_x);
		
		//interpolate fourth line
		arr_x[0] = f[x1][y4];
		arr_x[1] = f[x2][y4];
		arr_x[2] = f[x3][y4];
		arr_x[3] = f[x4][y4];
		arr_y[3] = cubic(xpos%1, arr_x);
		
		return cubic(ypos%1, arr_y);
	}
	
	
	
	public static float per_BiCub(float xpos, float ypos, float[][] f){
		int x = (int)xpos;
		int y = (int)ypos;
		
		int x1 = x-1;
		int x2 = x;
		int x3 = x+1;
		int x4 = x+2;
		
		int y1 = y-1;
		int y2 = y;
		int y3 = y+1;
		int y4 = y+2;
		
		/*
		if(x<1){x1=0;}
		if(y<1){y1=0;}
		
		if(x<0){x2=0; x3=1; x4=2;}
		if(y<0){y2=0; y3=1; y4=2;}
		
		if(x>ssx-1){x4=ssx+1;}
		if(y>ssy-1){y4=ssy+1;}
		
		if(x>ssx){x3=ssx+1;}
		if(y>ssy){y3=ssy+1;}
		
		if(x>ssx+1){x2=ssx+1; x1=ssx;}
		if(y>ssy+1){y2=ssy+1; y1=ssy;}
		*/
		
		if(y<1){y1=0;}
		if(y<0){y2=0; y3=1; y4=2;}
		if(y>ssy-1){y4=ssy+1;}
		if(y>ssy){y3=ssy+1;}
		if(y>ssy+1){y2=ssy+1; y1=ssy;}
		
		
		// Periodic X-Boarder
		if(x1<1){					x1 = ssx + x1;	}
		if(x2<1){					x2 = ssx + x2;	}
		if(x3<1){					x3 = ssx + x3;	}
		if(x4<1){					x4 = ssx + x4;	}
		if(x1>ssx){					x1 = (x1-ssx);	}			
		if(x2>ssx){					x2 = (x2-ssx);	}
		if(x3>ssx){					x3 = (x3-ssx);	}			
		if(x4>ssx){					x4 = (x4-ssx);	}
		
		
		
		
		float[] arr_x = new float[4];
		float[] arr_y = new float[4];
		
		//interpolate first line
		arr_x[0] = f[x1][y1];
		arr_x[1] = f[x2][y1];
		arr_x[2] = f[x3][y1];
		arr_x[3] = f[x4][y1];
		arr_y[0] = cubic(xpos%1, arr_x);
		
		//interpolate second line
		arr_x[0] = f[x1][y2];
		arr_x[1] = f[x2][y2];
		arr_x[2] = f[x3][y2];
		arr_x[3] = f[x4][y2];
		arr_y[1] = cubic(xpos%1, arr_x);
		
		//interpolate third line
		arr_x[0] = f[x1][y3];
		arr_x[1] = f[x2][y3];
		arr_x[2] = f[x3][y3];
		arr_x[3] = f[x4][y3];
		arr_y[2] = cubic(xpos%1, arr_x);
		
		//interpolate fourth line
		arr_x[0] = f[x1][y4];
		arr_x[1] = f[x2][y4];
		arr_x[2] = f[x3][y4];
		arr_x[3] = f[x4][y4];
		arr_y[3] = cubic(xpos%1, arr_x);
		
		return cubic(ypos%1, arr_y);
	}
	
	
	public static float per_BiLin(float xpos, float ypos, float[][] f){
		x1 = (int) xpos;
		x2 = x1+1;
		y1 = (int) ypos;
		y2 = y1+1;
		
		// Periodic X-Boarder
		if(x1<1){					x1 = ssx + x1;	}
		if(x2<1){					x2 = ssx + x2;	}
		if(x1>ssx){					x1 = (x1-ssx);	}			
		if(x2>ssx){					x2 = (x2-ssx);	}
		
		// Y-Boarder
		if(y1<0){					y1 = 0;	}
		if(y2<0){					y2 = 0;}
		if(y1>ssy+1){				y1 = ssy+1;}			
		if(y2>ssy+1){				y2 = ssy+1;}
		
		// Distances
		dx = (float) (xpos-Math.floor(xpos));//xpos%1;
		dy = (float) (ypos-Math.floor(ypos));//ypos%1;
		
		return  	(1-dx) *( 	f[x1][y1]*(1-dy) 	+ 	f[x1][y2]*dy  	)
					+  dx  *(	f[x2][y1]*(1-dy)	+  	f[x2][y2]*dy  	);
	}
	
	
	

	

	
	
}
