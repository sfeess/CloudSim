import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;

import javax.swing.*;


public class FluidViewer implements ActionListener, MouseListener,MouseMotionListener {

	private static final long serialVersionUID = 1L;
	
	static BufferedImage img;
	static Graphics2D onimg;
	static float[] pixelField,u,v;
	
	static FluidSolver fs = new FluidSolver();
	private static FluidPanel fp;
	
	static int ssx, ssy, sx, sy;
	static float scaleOut,dt;
	
	//Menu Stuff
	JMenuItem reset;
	JMenuItem close;
	JMenuItem values;
	JMenuItem vectors;
	JMenuItem stepcount;
	JMenuItem vort;
	JMenuItem paintVel,paintDens;
	static JSlider dtSlider ;
	
	static boolean dispVec,dispVal,dispVort,dispSteps,mDens,mVel;
	
	JPanel contentPanel;
	
	static int mx, my, mxOld, myOld;

	
	public static void init(){

		mx=my=myOld=mxOld=0;
		
		ssy = 100;
		ssx = 100;
		scaleOut = 3;
		
		sx = (int)scaleOut*ssx;
		sy = (int)scaleOut*ssy;
		
		//setup FluidSolver size
		fs.setup(ssx, ssy, 0.19F, scaleOut);
		fp = new FluidPanel(fs);
		
		//setup Output
		pixelField = new float[sx*sy];
		img = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		onimg = img.createGraphics();
		dispVec=mVel=true;
		dispVal=dispSteps=mDens=dispVort=false;
	}
	
	// MAIN METHOD
	public static void main(String[] args) {
		init();
		if(fs.step==0){
			try{
				//avoid refresh flicker
                Thread.sleep(10); //40= 25frames/1sec * 1000 millisec/sec
            }
            catch (InterruptedException e){}
		}
			
		new FluidViewer();
		while(fs.step<200000){
			fs.dt=dtSlider.getValue()/100f;

			//fs.vOld=Field.constField(ssx, ssy, (float) 0.01);
			//fs.uOld=Field.boxField(ssx, ssy, (float) Math.sin(fs.step/10));
			
			fs.step();
			
			
			
			fp.repaint();
			
		}
	}	
	
	// Constructor
	public FluidViewer(){
		JFrame frame = new JFrame("Fluid Viewer");
		JPanel topPanel = new JPanel();
		JPanel centerPanel = new JPanel();
		JPanel bottomPanel = new JPanel();
		
		frame.setLocation(100,100);
		frame.setSize(fs.sx+20,fs.sy+100);
		frame.setResizable(false);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		
		// Top Panel
		JMenuBar menuBar = new JMenuBar();
        JMenu file = new JMenu("File");
	    JMenu display = new JMenu("Display");
	    JMenu action = new JMenu("Action");
        
	    reset = new JMenuItem("Reset");
        reset.addActionListener(this);
        close = new JMenuItem("Close");
        close.addActionListener(this);
        
        vectors = new JMenuItem("Vectors");
        vectors.addActionListener(this);
        values = new JMenuItem("Velocity Values");
        values.addActionListener(this);
        stepcount = new JMenuItem("Frames");
        stepcount.addActionListener(this);
        vort = new JMenuItem("Vorticity");
        vort.addActionListener(this);
        
        paintVel = new JMenuItem("paint Velocity");
        paintVel.addActionListener(this);
        paintDens = new JMenuItem("paint Density");
        paintDens.addActionListener(this);
        
        menuBar.add(file);
        menuBar.add(action);
        menuBar.add(display);
        
        file.add(reset);
        file.add(close);
        display.add(vectors);
        display.add(values);
        display.add(stepcount);
        display.add(vort);
        action.add(paintVel);
        action.add(paintDens);
        
        topPanel.add(menuBar, BorderLayout.WEST);
         
        //CenterPanel 
       // FluidPanel fp = new FluidPanel();
        //centerPanel.setBackground(Color.blue);
        //centerPanel.add(fp);
        
        
         
		// Bottom Control Panel
        dtSlider = new JSlider(0,100);
        JLabel dtLabel = new JLabel("Timestep 0-1: ");
         
         bottomPanel.add(dtLabel);
         bottomPanel.add(dtSlider);
         
        fp.addMouseListener(this);
        fp.addMouseMotionListener(this);
        //centerPanel.setOpaque(false);  
        //frame.setGlassPane(centerPanel);
        //centerPanel.setVisible(true);  

        
        
        //frame.add(new FluidPanel(), BorderLayout.CENTER); 
        frame.add(topPanel, BorderLayout.NORTH);
        frame.add(fp, BorderLayout.CENTER);
        frame.add(bottomPanel, BorderLayout.SOUTH);
        
        
       
        
        frame.setVisible(true);
	}

	
	
	// linearisierung NUR FÜR PIXELFIELD!!! nicht sim grid
	public static int plin(int i, int j){
		return ((i)+(sx)*(j));
	}
	public static int flin(int i, int j){
		return ((i)+(ssx+2)*(j));
	}

	@Override
	public void actionPerformed(ActionEvent object) {
		// TODO Auto-generated method stub
		if (object.getSource() == reset){
			fs.reset();
       }
       if (object.getSource() == close){
    	   System.exit(0);
       }
       if (object.getSource() == values){
           dispVal=!dispVal; 
       }
       if (object.getSource() == vectors){
    	   dispVec=!dispVec; 
       }
       if (object.getSource() == stepcount){
    	   dispSteps=!dispSteps; 
       }
       if (object.getSource() == vort){
    	   dispVort=!dispVort;
       }
       if (object.getSource() == paintVel){
    	   mVel=!mVel;
       }
       if (object.getSource() == paintDens){
    	   mDens=!mDens;
       }
       
     
	}

	static boolean mouseDown = false;
	
	@Override
	public void mouseClicked(MouseEvent obj) {
		// TODO Auto-generated method stub
		
		
	}

	@Override
	public void mouseEntered(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void mouseExited(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void mousePressed(MouseEvent obj) {
		// TODO Auto-generated method stub

		mx=my=0;
		mxOld = (int) (obj.getX()/scaleOut);
		myOld = (int) ((sy-obj.getY())/scaleOut);
		
		
		
	}

	@Override
	public void mouseReleased(MouseEvent arg0) {
		// TODO Auto-generated method stub
		mouseDown=false;
		
	}

	@Override
	public void mouseDragged(MouseEvent obj) {
		//System.out.println("mouse dragged");
		mx = (int) (obj.getX()/scaleOut);
		my = (int) ((sy-obj.getY())/scaleOut);
		
		//System.out.println("mx "+mx);
		//System.out.println("my "+my);
		
		if(mDens && mx>2 && my>2 && ssx-2>mx && (ssy-2)>my){
			fs.dOld[flin(mx,my)]=1;
			fs.dOld[flin(mx-1,my)]=1;
			fs.dOld[flin(mx+1,my)]=1;
			fs.dOld[flin(mx,my-1)]=1;
			fs.dOld[flin(mx,my+1)]=1;
			fs.dOld[flin(mx-2,my)]=1;
			fs.dOld[flin(mx+2,my)]=1;
			fs.dOld[flin(mx,my-2)]=1;
			fs.dOld[flin(mx,my+2)]=1;
			
			fs.dOld[flin(mx-1,my+1)]=0.5f;
			fs.dOld[flin(mx+1,my+1)]=0.5f;
			fs.dOld[flin(mx-1,my-1)]=0.5f;
			fs.dOld[flin(mx+1,my-1)]=0.5f;
		}
		
		if(mVel && mx>2 && my>2 && ssx-2>mx && (ssy-2)>my){
			
		fs.uOld[flin(mx,my)]=mx-mxOld;
		fs.vOld[flin(mx,my)]=my-myOld;
		
		fs.uOld[flin(mx+1,my)]=mx-mxOld;
		fs.uOld[flin(mx,my+1)]=mx-mxOld;
		fs.uOld[flin(mx-1,my)]=mx-mxOld;
		fs.uOld[flin(mx,my-1)]=mx-mxOld;
		
		fs.vOld[flin(mx+1,my)]=my-myOld;
		fs.vOld[flin(mx,my+1)]=my-myOld;
		fs.vOld[flin(mx-1,my)]=my-myOld;
		fs.vOld[flin(mx,my-1)]=my-myOld;
		
		myOld=my;
		mxOld=mx;
		}
		
	}

	@Override
	public void mouseMoved(MouseEvent obj) {
		// TODO Auto-generated method stub
		
		
		
	}
}
		