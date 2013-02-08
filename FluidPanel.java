import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;

import javax.swing.*;


public class FluidPanel extends JPanel{

	private static final long serialVersionUID = 1L;
	
	static BufferedImage img;
	static Graphics2D onimg;
	static float[][] pixelField1,pixelField2, curlField, u,v;
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
		
		curlField = new float[sx][sy];
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
					pixelField1 = evaluate(FluidViewer.fs.qv);
					mapScale= 1f;
					pixelField2 = evaluate(FluidViewer.fs.qc);
					
				if(FluidViewer.dispVort)
					curlField=evaluate(FluidViewer.fs.vorticity);
				
				onimg.setColor(Color.red);
				u = evaluate(FluidViewer.fs.u);
				v = evaluate(FluidViewer.fs.v);
				
				WritableRaster raster= img.getRaster();
				
				// output Pixel Field
				for(int i=0; i<sx; i++){
					for(int j=0; j<sy; j++){
					
						int v =(int)((5*255*curlField[i][sy-1-j]));
						v = v>254 ? 255:v; 
						
						// cloud vapor
						qv =(int)((4*pixelField1[i][sy-1-j])*255);
						qv = qv<0 ? 0:qv; 
						
						//cloud water
						qc =(int)((110*pixelField2[i][sy-1-j])*255);
						qc = qc<0 ? 0:qc; 
						
					
						
						// Background Color
					 	
						int[] d = new int[3];	
						d[0]=51; d[1]=102; d[2]=153; 
						
						//qc=0;
						
						
						//d[0]=d[1]=d[2]=0; 
						
						d[0] = Math.min(255, d[0]+qc);
						d[1] = Math.min(255, d[1]+qc+qv);
						d[2] = Math.min(255, d[2]+qc);
						
						
						if(FluidViewer.dispVort){	d[0]=d[1]=d[2]=v; }	
				
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
	
	
	

}
		