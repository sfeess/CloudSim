import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;

import javax.swing.*;


public class FluidViewer extends JPanel {

	private static final long serialVersionUID = 1L;
	
	static BufferedImage img;
	static Graphics2D onimg;
	static float[] pixelField;
	static float[] vectorField;
	static float[] u,v;
	
	static FluidSolver fs = new FluidSolver();
	
	static int t=0;
	
	// Output PixelSize
	static int sx, sy;
	// Simulation grid size
	static int ssx, ssy;
	static int flip;
	
	// linearisierung NUR FÜR PIXELFIELD!!! nicht sim grid
	public static int plin(int i, int j){
		return ((i)+(sx)*(j));
	}
	public static int flin(int i, int j){
		return ((i)+(ssx+2)*(j));
	}
	
	public static void init(){

		
		sy = 200;
		sx = 200;
		ssx = 200;
		ssy = 200;
		
		//setup FluidSolver size
		fs.setup(ssx, ssy, 0.5F, sx, sy);
		
		
		//setup Output Field
		pixelField = new float[sx*sy];
		
		// setup BufferedImage
		img = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		onimg = img.createGraphics();
	}
	
	 public void update(Graphics g){ paint(g); }
	
	// MAIN METHOD
	public static void main(String[] args) {
		init();
		FluidViewer fw = new FluidViewer(sx,sy);
		
		while(true){

			fs.step();
			try
            {
                Thread.sleep(30);
            }
            catch (InterruptedException e)
            {
            }
            
			fw.repaint();
				//System.out.println("Size des arrays"+fs.out.length);
				//System.out.println("len 62 62 = "+lin(62,62,60));
				
		}
	}
	
	// Constructor
	public FluidViewer(int sizex, int sizey){
		JFrame frame = new JFrame("Fluid Viewer");
		frame.setLocation(100,100);
		frame.setSize(sizex+50,sizey+50);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.add(this);
		frame.setVisible(true);
	}

	// Paint Panel
	@Override
	public void paintComponent(Graphics g) {
		super.paintComponent(g);
	
		// BiLinear interpolation SimGrid to Pixels
		pixelField=fs.evaluate(sx,sy,fs.out);
		
		onimg.setColor(Color.red);
		u = fs.evaluate(sx,sy,fs.u);
		v = fs.evaluate(sx,sy,fs.v);
		
		// output Pixel Field
		WritableRaster raster= img.getRaster();
		
		for(int i=0; i<sx; i++){
			for(int j=0; j<sy; j++){
				int c =(int)((1-pixelField[plin(i,j)])*255);
				//c=(int)((((1-v[plin(i,j)])+1)/2)*255);
				int[] d ={c,c,c};
				raster.setPixel(i,j,d);	
			}
		}
		
		// Data for Vectors
		
		
		for(int i=5; i<sx; i+=10){
			for(int j=5; j<sy; j+=10){
				int u1 = (int) (7 * u[plin(i,j)] );
				int v1 = (int) (7* v[plin(i,j)]);
				// positive x -> to right**********positive y -> up
				if(i%50==5&&j%50==5){
					onimg.setColor(Color.gray);
					onimg.drawString("u="+(int)(100*u[plin(i,j)]),i,j);
					onimg.drawString("v="+(int)(100*v[plin(i,j)]),i,j+10);
				}
				onimg.setColor(Color.red);
				onimg.drawLine(i, j, i+u1, j+v1);
				
			}	
		}
		
		
		
		
		// image on Panel
		g.drawImage(img,0,0,this);
	}

}
		