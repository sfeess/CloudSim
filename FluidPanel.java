import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;

import javax.swing.*;


public class FluidPanel extends JPanel{

	private static final long serialVersionUID = 1L;
	
	static BufferedImage img;
	static Graphics2D onimg;
	static float[] pixelField,curlField, u,v;
	static FluidSolver fs;
	static int sx,sy,ssx,ssy;
	
	public FluidPanel(FluidSolver f){
		//fs=f;
		sx=FluidViewer.sx;
		sy=FluidViewer.sy;
		ssx=FluidViewer.ssx;
		ssy=FluidViewer.ssy;
		
		
		pixelField = new float[sx*sy];
		curlField = new float[sx*sy];
		img = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		onimg = img.createGraphics();

	}
	
	

	public void update(Graphics g){ paint(g); }

	// Paint Panel
	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
		
	
		// BiLinear interpolation SimGrid to Pixels
		pixelField=FluidViewer.fs.evaluate(sx,sy,FluidViewer.fs.d);
		curlField=FluidViewer.fs.evaluate(sx,sy,FluidViewer.fs.vorticity);
		
		onimg.setColor(Color.red);
		u = FluidViewer.fs.evaluate(sx,sy,FluidViewer.fs.u);
		v = FluidViewer.fs.evaluate(sx,sy,FluidViewer.fs.v);
		
		WritableRaster raster= img.getRaster();
		
		// output Pixel Field
		
		for(int i=0; i<sx; i++){
			for(int j=0; j<sy; j++){
			
				int v =(int)((10*255*curlField[plin(i,sy-1-j)]));
				v = v>254 ? 255:v; 
				
				//invert Y output
				int c =(int)((1-pixelField[plin(i,sy-1-j)])*255);
				c = c<0 ? 0:c; 
				
				int[] d ={c,c,c};
				
				if(FluidViewer.dispVort){	d[0]=d[1]=d[2]=v; }	
				
				raster.setPixel(i,j,d);	
			}
		}
		
		// paint Vectors
		int u1,v1;
			
		for(int i=5; i<sx; i+=10){
			for(int j=5; j<sy; j+=10){
					
				if(FluidViewer.dispVec){
				u1 = (int) (7 * u[plin(i,sy-1-j)] );
				v1 = (int) (-7* v[plin(i,sy-1-j)]);
				onimg.setColor(Color.red);
				onimg.drawLine(i, j, i+u1, j+v1);
				}
					
				if(i%100==35&&j%100==35&&FluidViewer.dispVal){
					onimg.setColor(Color.gray);
					onimg.drawString("u="+(u[plin(i,j)]),i,j);
					onimg.drawString("v="+(v[plin(i,j)]),i,j+10);
				}
			}	
		}
		
		
		if(FluidViewer.dispSteps){
		// Display frame number
		onimg.setColor(Color.black);
		onimg.drawString("Frame:"+FluidViewer.fs.step,8,15);
		}
		// image on Panel
		g.drawImage(img,0,0,this);
	}

	
	// linearisierung NUR FÜR PIXELFIELD!!! nicht sim grid
	public static int plin(int i, int j){
		return ((i)+(sx)*(j));
	}
	public static int flin(int i, int j){
		return ((i)+(ssx+2)*(j));
	}

}
		