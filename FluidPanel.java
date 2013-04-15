import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;

import javax.swing.*;


public class FluidPanel extends JPanel{

	private static final long serialVersionUID = 1L;
	
	static BufferedImage img;
	static Graphics2D onimg;
	static FluidSolver fs;
	static float scaleOut;
	static float[][] out_1,out_2, out_u, out_v, out_pt, out_d;
	static int c_qc, c_qv,c_pt,c_u,c_v, c_vort, c_d;
	FluidSolver f;
	
	
	public FluidPanel(FluidSolver f_solver, float scale){
		f=f_solver;
		scaleOut= scale;
		
		out_1 = new float[f.sx][f.sy];
		out_2 = new float[f.sx][f.sy];
		out_u = new float[f.sx][f.sy];
		out_v = new float[f.sx][f.sy];
		out_pt = new float[f.sx][f.sy];
		out_d = new float[f.sx][f.sy];
		
		img = new BufferedImage(f.sx, f.sy, BufferedImage.TYPE_INT_RGB);
		onimg = img.createGraphics();
		

	}
	
	

	public void update(Graphics g){ 
		paint(g); }

	
	public void outFields(){
		
		if(FluidViewer.dispMain==0 ){	
			out_1  = evaluate(f.qc);
			out_2  = evaluate(f.qv);
		}	
		else if(FluidViewer.dispMain==1){
			out_pt = evaluate(f.pt);
		}
		else if(FluidViewer.dispMain==2){
			out_u  = evaluate(f.u);
			out_v  = evaluate(f.v);
		}
		else if(FluidViewer.dispMain==3){
			out_d = evaluate(f.d);
		}
		if(FluidViewer.dispVec){
			out_u  = evaluate(f.u);
			out_v  = evaluate(f.v);
		}
		

		WritableRaster raster= img.getRaster();
		int[] d = new int[3];
		
		for(int i=0; i<f.sx; i++){
			for(int j=0; j<f.sy; j++){
			
				// Background Color
				//*********************************************************
			 	d[0]=51; d[1]=102; d[2]=153; 
				
				// Cloud Out
				//*********************************************************
				if(FluidViewer.dispMain==0){	
					c_qc =(int)((200*out_1[i][f.sy-1-j])*255);
					c_qc = c_qc<0 ? 0:c_qc; 
					
					if(FluidViewer.dispVapor){
						c_qv =(int)(20f*(out_2[i][f.sy-1-j])*255);
						c_qv = c_qv<0 ? 0:c_qv; 
					}
					else	c_qv=0;
					
					
					d[0] = Math.min(255, d[0]+c_qc);
					d[1] = Math.min(255, d[1]+c_qc+c_qv);
					d[2] = Math.min(255, d[2]+c_qc);
				}
				
				// Temperature Out
				//*********************************************************
				if(FluidViewer.dispMain==1){
					c_pt = (int) out_pt[i][f.sy-1-j];
					c_pt = c_pt>254 ? 255:c_pt; 
					
					d[0]=heatColor(c_pt).getRed();
					d[1]=heatColor(c_pt).getGreen();
					d[2]=heatColor(c_pt).getBlue();
				}
				
				// Velocity Out
				//*********************************************************
				if(FluidViewer.dispMain==2){
					c_u =(int)(1f*(out_u[i][f.sy-1-j])*255);
					c_u = Math.abs(c_u);
					c_u = c_u<0 ? 0:c_u; 
					
					c_v =(int)(1f*(out_v[i][f.sy-1-j])*255);
					c_v = Math.abs(c_v);
					c_v = c_v<0 ? 0:c_v; 
					
					
					d[0] = Math.min(255, c_u);
					d[1] = Math.min(255, c_v);
					d[2] = Math.min(255, 0);
				}
				
				// Density Out
				//*********************************************************
				if(FluidViewer.dispMain==3){
					c_d =(int)(1f*(out_d[i][f.sy-1-j])*255);
					c_d = c_d<0 ? 0:c_d; 
					
					
					d[0] = Math.min(255, c_d);
					d[1] = Math.min(255, c_d);
					d[2] = Math.min(255, c_d);
				}
				
				raster.setPixel(i,j,d);			
			}
		}
		
		
		
	}
	// Paint Panel
	@Override
	public void paintComponent(Graphics g) {
		
		super.paintComponent(g);
		
		

		// paint Vectors
		
		if(FluidViewer.dispVec){
			int u1,v1;
			
			for(int i=5; i<f.sx; i+=10){
				for(int j=5; j<f.sy; j+=10){
					
					u1 = (int) (out_u[i][f.sy-1-j]*7);
					v1 = (int) (out_v[i][f.sy-1-j]*-7);
					
					onimg.setColor(Color.red);
					onimg.drawLine(i, j, i+u1, j+v1);
					
					}
				}	
		}
		
		
		if(FluidViewer.dispSteps){
		// Display frame number
		onimg.setColor(Color.darkGray);
		onimg.drawString("Frame:"+f.step,9,16);
		onimg.setColor(Color.white);
		onimg.drawString("Frame:"+f.step,8,15);
		}
		// image on Panel
		g.drawImage(img,0,0,this);
		
	}

	
	
	
	
	
	
	
	public float[][] evaluate(float[][] field){
		
		//Interpolation.setSpace(f.sx, f.sy);
		float[][] mapped = new float[f.sx][f.sy];
		
		for(int i=0; i<f.sx; i++){
			for(int j=0; j<f.sy; j++){
				//mapped[i][j]= interpolate((float)(i/scaleOut+1),(float) (j/scaleOut+1),field); 
				mapped[i][j] = Interpolation.biCubic((float)(i/scaleOut+1), (float) (j/scaleOut+1), field); 
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
		