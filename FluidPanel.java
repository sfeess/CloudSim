import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;

import javax.swing.*;


public class FluidPanel extends JPanel{

	private static final long serialVersionUID = 1L;
	
	static BufferedImage img;
	static Graphics2D onimg;
	static float[][] pixelField1,pixelField2, ptField, u,v;
	static FluidSolver fs;
	static int sx,sy,ssx,ssy;
	static float scaleOut;
	static float mapScale;
	static int qc,qv;
	
	public FluidPanel(FluidSolver f){
		//fs=f;
		sx=FluidViewer.sx;
		sy=FluidViewer.sy;
		ssx=FluidViewer.ssx;
		ssy=FluidViewer.ssy;
		
		scaleOut= FluidViewer.scaleOut;
		
		
		pixelField1 = new float[sx][sy];
		pixelField2 = new float[sx][sy];
		
		ptField = new float[sx][sy];
		img = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		onimg = img.createGraphics();
		

	}
	
	

	public void update(Graphics g){ paint(g); }

	// Paint Panel
	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		
	
		// BiLinear interpolation SimGrid to Pixels
				if(!FluidViewer.dispVort)
					//pixelField1 = evaluate(FluidViewer.fs.qv);
					pixelField1 = evaluate(FluidViewer.fs.d);
					mapScale= 1f;
					//pixelField2 = evaluate(FluidViewer.fs.qc);
					
				if(FluidViewer.dispVort)
					ptField=evaluate(FluidViewer.fs.pt);
				
				onimg.setColor(Color.red);
				u = evaluate(FluidViewer.fs.u);
				v = evaluate(FluidViewer.fs.v);
				
				WritableRaster raster= img.getRaster();
				
				// output Pixel Field
				for(int i=0; i<sx; i++){
					for(int j=0; j<sy; j++){
					
						int v = (int) ptField[i][sy-1-j];
						//v = v>254 ? 255:v; 
						
						// cloud vapor
						qv =(int)((1f*pixelField1[i][sy-1-j])*255);
						qv = qv<0 ? 0:qv; 
						
						//cloud water
						//qc =(int)((100*pixelField2[i][sy-1-j])*255);
						//qc = qc<0 ? 0:qc; 
						
					
						
						// Background Color
					 	
						int[] d = new int[3];	
						d[0]=51; d[1]=102; d[2]=153; 
						//qc=0;
						
						
						//d[0]=d[1]=d[2]=0; 
						
						d[0] = Math.min(255, d[0]+qc);
						d[1] = Math.min(255, d[1]+qc+qv);
						d[2] = Math.min(255, d[2]+qc);
						
						d[0] = Math.min(255, qv);
						d[1] = Math.min(255, qv);
						d[2] = Math.min(255, qv);
						
						if(FluidViewer.dispVort){	
							d[0]=heatColor(v).getRed();
							d[1]=heatColor(v).getGreen();
							d[2]=heatColor(v).getBlue();
						}	
				
						raster.setPixel(i,j,d);	
			}
		}
		
		// paint Vectors
		int u1,v1;
			
		for(int i=5; i<sx; i+=10){
			for(int j=5; j<sy; j+=10){
					
				if(FluidViewer.dispVec){
				u1 = (int) (7 * u[i][sy-1-j] );
				v1 = (int) (-7* v[i][sy-1-j]);
				onimg.setColor(Color.red);
				onimg.drawLine(i, j, i+u1, j+v1);
				}
					
				if(i%100==35&&j%100==35&&FluidViewer.dispVal){
					onimg.setColor(Color.gray);
					onimg.drawString("u="+(u[i][j]),i,j);
					onimg.drawString("v="+(v[i][j]),i,j+10);
				}
			}	
		}
		
		
		if(FluidViewer.dispSteps){
		// Display frame number
		onimg.setColor(Color.magenta);
		onimg.drawString("Frame:"+FluidViewer.fs.step,8,15);
		}
		// image on Panel
		g.drawImage(img,0,0,this);
		
	}

	
	
	
	
	
	
	
	public float[][] evaluate(float[][] f){
		float[][] mapped = new float[sx][sy];
		
		for(int i=0; i<sx; i++){
			for(int j=0; j<sy; j++){
				mapped[i][j]= interpolate((float)(i/scaleOut+1),(float) (j/scaleOut+1),f); //f[(int) (i/scaleOut+1)][(int) (j/scaleOut+1)];	
			}
		}
		return mapped;
	}
	
	
	public float interpolate(float xpos, float ypos, float[][] f){
		// Sample positions
		int x1 = (int) Math.floor(xpos);
		int x2 = x1+1;
		int y1 = (int) Math.floor(ypos);
		int y2 = y1+1;
		// Distances
		float dx = (float) (xpos-Math.floor(xpos));
		float dy = (float) (ypos-Math.floor(ypos));
		// Interpolated Value
		return  	(1-dx) *( 	f[x1][y1]*(1-dy) 	+ 	f[x1][y2]*dy  	)
					+  dx  *(	f[x2][y1]*(1-dy)	+  	f[x2][y2]*dy  	);
		

		
	}
	
	public static Color heatColor(float temp){
		
		Color c1 = new Color(50, 0, 255);
		Color c2 = new Color(0, 230, 160);
		Color c3 = new Color(255, 255, 0);
		Color c4 = new Color(180, 0, 50);
		Color c5 = new Color(255, 0, 255);
		
		float  p[] = {200f,250f,300f,350f,450};
		
		if(temp<p[1]){
			c1 =  new Color(
							(int) (c1.getRed()    +  ( ((temp-p[0])/(p[1]-p[0]))  * (c2.getRed()  -c1.getRed() ) )), 
							(int) (c1.getGreen()  +  ( ((temp-p[0])/(p[1]-p[0]))  * (c2.getGreen()-c1.getGreen() ) )), 
							(int) (c1.getBlue()   +  ( ((temp-p[0])/(p[1]-p[0]))  * (c2.getBlue() -c1.getBlue() ) ))
							); 
		}
		else if(temp<p[2]){
			c1 =  new Color(
							(int) (c2.getRed()    +  ( ((temp-p[1])/(p[2]-p[1]))  * (c3.getRed()  -c2.getRed() ) )), 
							(int) (c2.getGreen()  +  ( ((temp-p[1])/(p[2]-p[1]))  * (c3.getGreen()-c2.getGreen() ) )), 
							(int) (c2.getBlue()   +  ( ((temp-p[1])/(p[2]-p[1]))  * (c3.getBlue() -c2.getBlue() ) ))
							); 
		}
		else if(temp<p[3]){
			c1 =  new Color(
							(int) (c3.getRed()    +  ( ((temp-p[2])/(p[3]-p[2]))  * (c4.getRed()-c3.getRed() ) )),
							(int) (c3.getGreen()  +  ( ((temp-p[2])/(p[3]-p[2]))  * (c4.getGreen()-c3.getGreen() ) )),
							(int) (c3.getBlue()   +  ( ((temp-p[2])/(p[3]-p[2]))  * (c4.getBlue() -c3.getBlue() ) ))
							); 
		}
		else if(temp<p[4]){
			c1 =  new Color(
							(int) (c4.getRed()    +  ( ((temp-p[3])/(p[4]-p[3]))  * (c5.getRed()  -c4.getRed() ) )), 
							(int) (c4.getGreen()  +  ( ((temp-p[3])/(p[4]-p[3]))  * (c5.getGreen()-c4.getGreen() ) )), 
							(int) (c4.getBlue()   +  ( ((temp-p[3])/(p[4]-p[3]))  * (c5.getBlue() -c4.getBlue() ) ))
							); 
		}
		
		else c1= new Color(255,0,0);
		
		return c1;
		
	}
	
	
	

}
		