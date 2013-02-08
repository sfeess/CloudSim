import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;

import javax.swing.*;
import javax.swing.border.Border;





public class FluidViewer implements ActionListener, MouseListener,MouseMotionListener {

	static BufferedImage img;
	static Graphics2D onimg;
	static float[][] pixelField,u,v;
	
	static FluidSolver fs = new FluidSolver();
	private static FluidPanel fp;
	
	static int ssx, ssy, sx, sy;
	static float scaleOut,dt, fps;
	static long time;
	
	//Menu Stuff
	JMenuItem reset;
	JMenuItem close;
	JMenuItem values;
	JMenuItem debug;
	JMenuItem vectors;
	JMenuItem stepcount;
	JMenuItem vort;
	JMenuItem paintVel,paintDens;
	 static JLabel lblDebugvalue1;
	 static JLabel lblDebugvalue2;
	 static JLabel lblDebugvalue3;
	 static JLabel lblDebugvalue4;
	 static JLabel lblDebugvalue5;
	 static JLabel lblDebugvalue6;
	 static JLabel lblDebugPos;
	 
	JButton btnPlay = new JButton();
	JButton btnPause = new JButton();
	JButton btnStop = new JButton();
	
	static boolean dispVec,dispVal,dispVort,dispSteps,mDens,mVel,dispDebug, dispSettings, play, stop;
	
	JPanel contentPanel;
	
	
	static int mx, my, mxOld, myOld;

	
	
	
	// MAIN METHOD
	public static void main(String[] args) {
		init();
		new FluidViewer();
		
		// Simulation Loop
		//*****************************************************************
		while(fs.step<20000){
			if(play&&!stop){
				fps= (int)(1/((System.nanoTime()-time)/1000000000f)*10);
				fps /= 10f;
				time = System.nanoTime();
				
				if(fs.step<30000)
					fs.step();
				try{//avoid refresh flicker
					Thread.sleep(0); //40= 25frames/1sec * 1000 millisec/sec
				}
				catch (InterruptedException e){}
				
				fp.repaint();
				// Output Data
				//WriteData.imgOut();
				//WriteData.boxValuesOut(10,15,0,5);
			}
			if(dispDebug)refreshDebug(mx,my);
			
			
		}
		// End Simulation Loop
		// ************************************************************
	}	
	
	
	public static void init(){

		mx=my=myOld=mxOld=0;
		
		ssy = 100;
		ssx = 200;
		scaleOut = 3;
		
		sx = (int)scaleOut*ssx;
		sy = (int)scaleOut*ssy;
		//setup FluidSolver size
		fs.setup(ssx, ssy, 0.19F, scaleOut);
		fp = new FluidPanel(fs);
		fp.setBounds(5, 26, sx, sy);
		
		//setup Output
		pixelField = new float[sx][sy];
		img = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		onimg = img.createGraphics();
		dispVec=mVel=dispSteps=stop=true;
		dispVal=mDens=dispVort=dispDebug=dispVec=dispSettings=play=false;
		
		time = System.currentTimeMillis();
	}
	
	
	// Constructor
	public FluidViewer(){
		JFrame frame = new JFrame("Fluid Viewer");
		frame.setLocation(100,100);
		frame.setSize(sx+180,sy+190); 
		frame.getContentPane().setLayout(null);
		frame.setResizable(false);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().setBackground(new Color(105, 105, 105));
		
		frame.setIconImage(new ImageIcon("cloud-icon.png").getImage());
		
		JPanel topPanel = new JPanel();
		topPanel.setBounds(0, 0, sx+500, 21);
		frame.getContentPane().add(topPanel);
		
		JPanel bottomPanel = new JPanel();
		bottomPanel.setBackground(Color.DARK_GRAY);
		bottomPanel.setBounds(5, sy+31, sx, 125);
		bottomPanel.setLayout(null);
		
		frame.getContentPane().add(bottomPanel);
		
		
		btnPlay.setBounds((int)(sx/2)-56, 5, 34, 34);
		btnPlay.setIcon(new ImageIcon("img\\play.png"));
		btnPlay.setBorderPainted(false);
		btnPlay.addActionListener(this);
		btnPlay.addMouseListener(this);
		bottomPanel.add(btnPlay);
		
		btnPause.setBounds((int)(sx/2)-17, 5, 34, 34);
		btnPause.setIcon(new ImageIcon("img\\pause.png"));
		btnPause.setBorderPainted(false);
		btnPause.addActionListener(this);
		btnPause.addMouseListener(this);
		bottomPanel.add(btnPause);
		
		btnStop.setIcon(new ImageIcon("img\\stop.png"));
		btnStop.setBounds((int)(sx/2)+22, 5, 34, 34);
		btnStop.setBorderPainted(false);
		btnStop.addActionListener(this);
		btnStop.addMouseListener(this);
		bottomPanel.add(btnStop);
		 
		
		fp.addMouseListener(this);
        fp.addMouseMotionListener(this);
       
       
        
        
        
        JMenuBar menuBar = new JMenuBar();
        JMenu file = new JMenu("File");
        JMenu display = new JMenu("Display");
        JMenu action = new JMenu("Action");
        	
	    reset = new JMenuItem("Reset");
	    reset.addActionListener(this);
	    close = new JMenuItem("Close");
	    close.addActionListener(this);
	    
	    debug = new JMenuItem("Display Settings");
	    debug.addActionListener(this);
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
	    topPanel.setLayout(new GridLayout(0, 1, 0, 0));
	    
	    menuBar.add(file);
	    menuBar.add(action);
	    menuBar.add(display);
	    
	    file.add(reset);
	    file.add(close);
	    display.add(debug);
	    display.add(vectors);
	    display.add(values);
	    display.add(stepcount);
	    display.add(vort);
	    action.add(paintVel);
	    action.add(paintDens);
	    topPanel.add(menuBar);
	    frame.getContentPane().add(fp);
        
	        
       
       settingsPanel = new JPanel();
	   settingsPanel.setBackground(Color.DARK_GRAY);
	   settingsPanel.setBounds(sx+10, 26, 160, sy);
	        frame.getContentPane().add(settingsPanel);
	        GridBagLayout gbl_settingsPanel = new GridBagLayout();
	        gbl_settingsPanel.columnWidths = new int[] {60, 0};
	        gbl_settingsPanel.rowHeights = new int[]{20, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	        gbl_settingsPanel.columnWeights = new double[]{1.0, 1.0};
	        gbl_settingsPanel.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Double.MIN_VALUE};
	        settingsPanel.setLayout(gbl_settingsPanel);
	        
	        JLabel lblInitialValuesrestart = new JLabel("Initial Values");
	        lblInitialValuesrestart.setForeground(Color.WHITE);
	        lblInitialValuesrestart.setFont(new Font("Monospaced", Font.PLAIN, 11));
	        GridBagConstraints gbc_lblInitialValuesrestart = new GridBagConstraints();
	        gbc_lblInitialValuesrestart.anchor = GridBagConstraints.WEST;
	        gbc_lblInitialValuesrestart.gridwidth = 2;
	        gbc_lblInitialValuesrestart.insets = new Insets(1, 5, 5, 0);
	        gbc_lblInitialValuesrestart.gridx = 0;
	        gbc_lblInitialValuesrestart.gridy = 0;
	        settingsPanel.add(lblInitialValuesrestart, gbc_lblInitialValuesrestart);
	        
	        txtAlt = new JTextField();
	        txtAlt.setForeground(Color.WHITE);
	        txtAlt.setBackground(Color.GRAY);
	        txtAlt.setHorizontalAlignment(SwingConstants.LEFT);
	        txtAlt.setText("6000");
	        GridBagConstraints gbc_txtAlt = new GridBagConstraints();
	        gbc_txtAlt.fill = GridBagConstraints.HORIZONTAL;
	        gbc_txtAlt.insets = new Insets(5, 5, 5, 5);
	        gbc_txtAlt.gridx = 0;
	        gbc_txtAlt.gridy = 1;
	        settingsPanel.add(txtAlt, gbc_txtAlt);
	        txtAlt.setColumns(10);
	        txtAlt.addActionListener(this);
	        
	        lblMaxAltitude = new JLabel("max. Altitude");
	        lblMaxAltitude.setFont(new Font("Monospaced", Font.PLAIN, 11));
	        lblMaxAltitude.setForeground(Color.LIGHT_GRAY);
	        GridBagConstraints gbc_lblMaxAltitude = new GridBagConstraints();
	        gbc_lblMaxAltitude.insets = new Insets(5, 0, 5, 0);
	        gbc_lblMaxAltitude.anchor = GridBagConstraints.WEST;
	        gbc_lblMaxAltitude.gridx = 1;
	        gbc_lblMaxAltitude.gridy = 1;
	        settingsPanel.add(lblMaxAltitude, gbc_lblMaxAltitude);
	        
	        txtHum = new JTextField();
	        txtHum.setForeground(Color.WHITE);
	        txtHum.setBackground(Color.GRAY);
	        txtHum.setText("0.5");
	        GridBagConstraints gbc_txtHum = new GridBagConstraints();
	        gbc_txtHum.fill = GridBagConstraints.HORIZONTAL;
	        gbc_txtHum.insets = new Insets(0, 5, 5, 5);
	        gbc_txtHum.gridx = 0;
	        gbc_txtHum.gridy = 2;
	        settingsPanel.add(txtHum, gbc_txtHum);
	        txtHum.setColumns(10);
	        txtHum.addActionListener(this);
	        
	        lblHumidity = new JLabel("Humidity");
	        lblHumidity.setFont(new Font("Monospaced", Font.PLAIN, 11));
	        lblHumidity.setForeground(Color.LIGHT_GRAY);
	        lblHumidity.setToolTipText("Defines the initial humidity of the air.");
	        GridBagConstraints gbc_lblHumidity = new GridBagConstraints();
	        gbc_lblHumidity.insets = new Insets(0, 0, 5, 0);
	        gbc_lblHumidity.anchor = GridBagConstraints.WEST;
	        gbc_lblHumidity.gridx = 1;
	        gbc_lblHumidity.gridy = 2;
	        settingsPanel.add(lblHumidity, gbc_lblHumidity);
	        
	        txtTlr = new JTextField();
	        txtTlr.setForeground(Color.WHITE);
	        txtTlr.setBackground(Color.GRAY);
	        txtTlr.setText("0.6");
	        GridBagConstraints gbc_txtTlr = new GridBagConstraints();
	        gbc_txtTlr.insets = new Insets(0, 5, 5, 5);
	        gbc_txtTlr.fill = GridBagConstraints.HORIZONTAL;
	        gbc_txtTlr.gridx = 0;
	        gbc_txtTlr.gridy = 3;
	        settingsPanel.add(txtTlr, gbc_txtTlr);
	        txtTlr.setColumns(10);
	        txtTlr.addActionListener(this);
	        
	        lblLapseRate = new JLabel("Lapse Rate");
	        lblLapseRate.setFont(new Font("Monospaced", Font.PLAIN, 11));
	        lblLapseRate.setForeground(Color.LIGHT_GRAY);
	        lblLapseRate.setToolTipText("Decrease of temperature per 100m");
	        GridBagConstraints gbc_lblLapseRate = new GridBagConstraints();
	        gbc_lblLapseRate.anchor = GridBagConstraints.WEST;
	        gbc_lblLapseRate.insets = new Insets(0, 0, 5, 0);
	        gbc_lblLapseRate.gridx = 1;
	        gbc_lblLapseRate.gridy = 3;
	        settingsPanel.add(lblLapseRate, gbc_lblLapseRate);
	        
	        lblInteractivValues = new JLabel("Interactive Values");
	        lblInteractivValues.setForeground(Color.WHITE);
	        lblInteractivValues.setFont(new Font("Monospaced", Font.PLAIN, 11));
	        GridBagConstraints gbc_lblInteractivValues = new GridBagConstraints();
	        gbc_lblInteractivValues.anchor = GridBagConstraints.WEST;
	        gbc_lblInteractivValues.gridwidth = 2;
	        gbc_lblInteractivValues.insets = new Insets(15, 5, 5, 0);
	        gbc_lblInteractivValues.gridx = 0;
	        gbc_lblInteractivValues.gridy = 4;
	        settingsPanel.add(lblInteractivValues, gbc_lblInteractivValues);
	        
	        txtVort = new JTextField();
	        txtVort.setForeground(Color.WHITE);
	        txtVort.setBackground(Color.GRAY);
	        txtVort.setHorizontalAlignment(SwingConstants.LEFT);
	        txtVort.setText("0.4");
	        GridBagConstraints gbc_txtVort = new GridBagConstraints();
	        gbc_txtVort.fill = GridBagConstraints.HORIZONTAL;
	        gbc_txtVort.insets = new Insets(0, 5, 5, 5);
	        gbc_txtVort.gridx = 0;
	        gbc_txtVort.gridy = 5;
	        settingsPanel.add(txtVort, gbc_txtVort);
	        txtVort.setColumns(10);
	        txtVort.addActionListener(this);
	        
	        lblVorticityConfinement = new JLabel("Vorticity");
	        lblVorticityConfinement.setFont(new Font("Monospaced", Font.PLAIN, 11));
	        lblVorticityConfinement.setForeground(Color.LIGHT_GRAY);
	        lblVorticityConfinement.setHorizontalAlignment(SwingConstants.CENTER);
	        lblVorticityConfinement.setToolTipText("controls the 'curlyness' of the simulation");
	        GridBagConstraints gbc_lblVorticityConfinement = new GridBagConstraints();
	        gbc_lblVorticityConfinement.insets = new Insets(0, 0, 5, 0);
	        gbc_lblVorticityConfinement.anchor = GridBagConstraints.WEST;
	        gbc_lblVorticityConfinement.gridx = 1;
	        gbc_lblVorticityConfinement.gridy = 5;
	        settingsPanel.add(lblVorticityConfinement, gbc_lblVorticityConfinement);
	        
	        txtBuoy = new JTextField();
	        txtBuoy.setForeground(Color.WHITE);
	        txtBuoy.setBackground(Color.GRAY);
	        txtBuoy.setHorizontalAlignment(SwingConstants.LEFT);
	        txtBuoy.setText("1.0");
	        GridBagConstraints gbc_txtBuoy = new GridBagConstraints();
	        gbc_txtBuoy.insets = new Insets(0, 5, 5, 5);
	        gbc_txtBuoy.fill = GridBagConstraints.HORIZONTAL;
	        gbc_txtBuoy.gridx = 0;
	        gbc_txtBuoy.gridy = 6;
	        settingsPanel.add(txtBuoy, gbc_txtBuoy);
	        txtBuoy.setColumns(10);
	        txtBuoy.addActionListener(this);
	        
	        lblBuoyancy = new JLabel("Buoyancy");
	        lblBuoyancy.setFont(new Font("Monospaced", Font.PLAIN, 11));
	        lblBuoyancy.setForeground(Color.LIGHT_GRAY);
	        lblBuoyancy.setToolTipText("Buoyancy defines the strength of warm air pushing up and cold air dragging down");
	        GridBagConstraints gbc_lblBuoyancy = new GridBagConstraints();
	        gbc_lblBuoyancy.insets = new Insets(0, 0, 5, 0);
	        gbc_lblBuoyancy.anchor = GridBagConstraints.WEST;
	        gbc_lblBuoyancy.gridx = 1;
	        gbc_lblBuoyancy.gridy = 6;
	        settingsPanel.add(lblBuoyancy, gbc_lblBuoyancy);
	        
	        txtDt = new JTextField();
	        txtDt.setForeground(Color.WHITE);
	        txtDt.setBackground(Color.GRAY);
	        txtDt.setText("0.5");
	        GridBagConstraints gbc_txtDt = new GridBagConstraints();
	        gbc_txtDt.insets = new Insets(0, 5, 5, 5);
	        gbc_txtDt.fill = GridBagConstraints.HORIZONTAL;
	        gbc_txtDt.gridx = 0;
	        gbc_txtDt.gridy = 7;
	        settingsPanel.add(txtDt, gbc_txtDt);
	        txtDt.setColumns(10);
	        txtDt.addActionListener(this);
	        
	        JLabel dtLabel = new JLabel("Timestep");
	        dtLabel.setToolTipText("the timestep used in one simulation frame");
	        GridBagConstraints gbc_dtLabel = new GridBagConstraints();
	        gbc_dtLabel.insets = new Insets(0, 0, 5, 0);
	        gbc_dtLabel.anchor = GridBagConstraints.WEST;
	        gbc_dtLabel.gridx = 1;
	        gbc_dtLabel.gridy = 7;
	        settingsPanel.add(dtLabel, gbc_dtLabel);
	        dtLabel.setFont(new Font("Monospaced", Font.PLAIN, 11));
	        dtLabel.setForeground(Color.LIGHT_GRAY);
			
			
			JPanel debugValues = new JPanel();
			debugValues.setBounds(sx+10, sy+31, 160, 125);
			frame.getContentPane().add(debugValues);
			debugValues.setBackground(Color.DARK_GRAY);
			GridBagLayout gbl_debugValues = new GridBagLayout();
			gbl_debugValues.columnWidths = new int[]{160, 0};
			gbl_debugValues.rowHeights = new int[]{17, 17, 17, 17, 17, 17, 17, 0};
			gbl_debugValues.columnWeights = new double[]{0.0, Double.MIN_VALUE};
			gbl_debugValues.rowWeights = new double[]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Double.MIN_VALUE};
			debugValues.setLayout(gbl_debugValues);
			lblDebugPos = new JLabel("x:000 y:000");
			lblDebugPos.setFont(new Font("Monospaced", Font.PLAIN, 11));
			lblDebugPos.setForeground(Color.LIGHT_GRAY);
			GridBagConstraints gbc_lblDebugPos = new GridBagConstraints();
			gbc_lblDebugPos.insets = new Insets(5, 5, 0, 0);
			gbc_lblDebugPos.fill = GridBagConstraints.BOTH;
			gbc_lblDebugPos.gridx = 0;
			gbc_lblDebugPos.gridy = 0;
			debugValues.add(lblDebugPos, gbc_lblDebugPos);
			lblDebugvalue1 = new JLabel("v1");
			lblDebugvalue1.setFont(new Font("Monospaced", Font.PLAIN, 11));
			lblDebugvalue1.setHorizontalAlignment(SwingConstants.LEFT);
			lblDebugvalue1.setForeground(Color.LIGHT_GRAY);
			GridBagConstraints gbc_lblDebugvalue1 = new GridBagConstraints();
			gbc_lblDebugvalue1.insets = new Insets(0, 5, 0, 0);
			gbc_lblDebugvalue1.fill = GridBagConstraints.BOTH;
			gbc_lblDebugvalue1.gridx = 0;
			gbc_lblDebugvalue1.gridy = 1;
			debugValues.add(lblDebugvalue1, gbc_lblDebugvalue1);
			lblDebugvalue2 = new JLabel("v2");
			lblDebugvalue2.setFont(new Font("Monospaced", Font.PLAIN, 11));
			lblDebugvalue2.setForeground(Color.LIGHT_GRAY);
			GridBagConstraints gbc_lblDebugvalue2 = new GridBagConstraints();
			gbc_lblDebugvalue2.insets = new Insets(0, 5, 0, 0);
			gbc_lblDebugvalue2.fill = GridBagConstraints.BOTH;
			gbc_lblDebugvalue2.gridx = 0;
			gbc_lblDebugvalue2.gridy = 2;
			debugValues.add(lblDebugvalue2, gbc_lblDebugvalue2);
			lblDebugvalue3 = new JLabel("v3");
			lblDebugvalue3.setFont(new Font("Monospaced", Font.PLAIN, 11));
			lblDebugvalue3.setForeground(Color.LIGHT_GRAY);
			GridBagConstraints gbc_lblDebugvalue3 = new GridBagConstraints();
			gbc_lblDebugvalue3.insets = new Insets(0, 5, 0, 0);
			gbc_lblDebugvalue3.fill = GridBagConstraints.BOTH;
			gbc_lblDebugvalue3.gridx = 0;
			gbc_lblDebugvalue3.gridy = 3;
			debugValues.add(lblDebugvalue3, gbc_lblDebugvalue3);
			lblDebugvalue4 = new JLabel("v4");
			lblDebugvalue4.setFont(new Font("Monospaced", Font.PLAIN, 11));
			lblDebugvalue4.setForeground(Color.LIGHT_GRAY);
			GridBagConstraints gbc_lblDebugvalue4 = new GridBagConstraints();
			gbc_lblDebugvalue4.insets = new Insets(0, 5, 0, 0);
			gbc_lblDebugvalue4.fill = GridBagConstraints.BOTH;
			gbc_lblDebugvalue4.gridx = 0;
			gbc_lblDebugvalue4.gridy = 4;
			debugValues.add(lblDebugvalue4, gbc_lblDebugvalue4);
			lblDebugvalue5 = new JLabel("v5");
			lblDebugvalue5.setFont(new Font("Monospaced", Font.PLAIN, 11));
			lblDebugvalue5.setForeground(Color.LIGHT_GRAY);
			GridBagConstraints gbc_lblDebugvalue5 = new GridBagConstraints();
			gbc_lblDebugvalue5.insets = new Insets(0, 5, 0, 0);
			gbc_lblDebugvalue5.fill = GridBagConstraints.BOTH;
			gbc_lblDebugvalue5.gridx = 0;
			gbc_lblDebugvalue5.gridy = 5;
			debugValues.add(lblDebugvalue5, gbc_lblDebugvalue5);
			lblDebugvalue6 = new JLabel("v6");
			lblDebugvalue6.setFont(new Font("Monospaced", Font.PLAIN, 11));
			lblDebugvalue6.setForeground(Color.LIGHT_GRAY);
			GridBagConstraints gbc_lblDebugvalue6 = new GridBagConstraints();
			gbc_lblDebugvalue6.insets = new Insets(0, 5, 0, 0);
			gbc_lblDebugvalue6.fill = GridBagConstraints.BOTH;
			gbc_lblDebugvalue6.gridx = 0;
			gbc_lblDebugvalue6.gridy = 6;
			debugValues.add(lblDebugvalue6, gbc_lblDebugvalue6);
		    
        
        
        
        
        frame.setVisible(true);
	}

	public static void refreshDebug(int x, int y){
		//x=40;		y=0;		
		lblDebugPos.setText("x: "+x+" y: "+y);
		lblDebugvalue1.setText("u:  "+FluidViewer.fs.u[x][y]);
		lblDebugvalue2.setText("v:  "+FluidViewer.fs.v[x][y]);
		lblDebugvalue3.setText("pt: "+FluidViewer.fs.pt[x][y]);
		lblDebugvalue4.setText("qc: "+FluidViewer.fs.qc[x][y]);
		lblDebugvalue5.setText("qv: "+FluidViewer.fs.qv[x][y]);
		lblDebugvalue6.setText("d:  "+FluidViewer.fs.d[x][y]);
	}
	
	// linearisierung NUR FÜR PIXELFIELD!!! nicht sim grid
	//public static int plin(int i, int j){	return ((i)+(sx)*(j));	}
	//public static int flin(int i, int j){	return ((i)+(ssx+2)*(j));	}

	@Override
	public void actionPerformed(ActionEvent object) {
		
		// Value Input
		//*****************************************************************
		
		if(object.getSource() == txtDt){
			fs.dt = textToFloat(txtDt.getText(), 0.001f, 1f, fs.dt, 1);
			txtDt.setText(String.valueOf(fs.dt));
		}
		else if(object.getSource() == txtTlr){
			fs.tlr = textToFloat(txtTlr.getText(), 0.55f, 0.99f, fs.tlr, 100);
			txtTlr.setText(String.valueOf(fs.tlr*100));
		}
		else if(object.getSource() == txtHum){
			fs.hum = textToFloat(txtHum.getText(), 0,1, fs.hum,1);
			txtHum.setText(String.valueOf(fs.hum));
		}
		else if(object.getSource() == txtBuoy){
			fs.buoyancy = textToFloat(txtBuoy.getText(), 0,5, fs.buoyancy,1);
			txtBuoy.setText(String.valueOf(fs.buoyancy));
		}
		else if(object.getSource() == txtVort){
			fs.vort = textToFloat(txtVort.getText(), 0,1, fs.vort,1);
			txtVort.setText(String.valueOf(fs.vort));
		}
		else if(object.getSource() == txtAlt){
			fs.maxAlt = textToFloat(txtAlt.getText(), 0,10000, fs.maxAlt,1);
			txtAlt.setText(String.valueOf(fs.maxAlt));
		}
		
		// Button Input
		//*****************************************************************
		else if(object.getSource() == btnPlay ){
			play = true;
			stop = false;
		}
		else if(object.getSource() == btnPause ){
			play = false;
			stop = false;
		}
		else if(object.getSource() == btnStop ){
			fs.reset();
			play = false;
			stop = true;
		}
		
		// Menu Input
		//*****************************************************************
		else if (object.getSource() == reset){
			fs.reset();
		}
		else if (object.getSource() == close){
    	   System.exit(0);
		}
		else if (object.getSource() == debug){
    	   dispDebug= !dispDebug;
          
		}
		else if (object.getSource() == values){
           dispVal=!dispVal; 
		}
		else if (object.getSource() == vectors){
    	   dispVec=!dispVec; 
		}
		else if (object.getSource() == stepcount){
    	   dispSteps=!dispSteps; 
		}
		else if (object.getSource() == vort){
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
	private JPanel settingsPanel;
	private JLabel lblVorticityConfinement;
	private JLabel lblBuoyancy;
	private JLabel lblHumidity;
	private JLabel lblMaxAltitude;
	private JTextField txtVort;
	private JTextField txtBuoy;
	private JTextField txtHum;
	private JTextField txtAlt;
	private JTextField txtTlr;
	private JLabel lblLapseRate;
	private JLabel lblInteractivValues;
	private JTextField txtDt;
	
	@Override
	public void mouseClicked(MouseEvent obj) {
		
		
	}

	@Override
	public void mouseEntered(MouseEvent object) {
		// Button animation
		//*************************************************
		if(object.getSource() == btnStop ){		btnStop.setIcon(new ImageIcon("img\\stop_hov.png"));	}
		if(object.getSource() == btnPlay ){		btnPlay.setIcon(new ImageIcon("img\\play_hov.png"));	}
		if(object.getSource() == btnPause ){	btnPause.setIcon(new ImageIcon("img\\pause_hov.png"));	}
		
	}

	@Override
	public void mouseExited(MouseEvent object) {
		// Button animation
		//*************************************************
		if(object.getSource() == btnStop ){		btnStop.setIcon(new ImageIcon("img\\stop.png"));	}
		if(object.getSource() == btnPlay ){		btnPlay.setIcon(new ImageIcon("img\\play.png"));	}
		if(object.getSource() == btnPause ){	btnPause.setIcon(new ImageIcon("img\\pause.png"));	}
		
	}

	@Override
	public void mousePressed(MouseEvent obj) {
		// Button animation
		//*************************************************
		if(obj.getSource() == btnStop ){		btnStop.setIcon(new ImageIcon("img\\stop_dwn.png"));	}
		if(obj.getSource() == btnPlay ){		btnPlay.setIcon(new ImageIcon("img\\play_dwn.png"));	}
		if(obj.getSource() == btnPause ){		btnPause.setIcon(new ImageIcon("img\\pause_dwn.png"));	}
		
		mx=my=0;
		mxOld = (int) (obj.getX()/scaleOut);
		myOld = (int) ((sy-obj.getY())/scaleOut);
		
		
		
	}

	@Override
	public void mouseReleased(MouseEvent object) {
		// Button animation
		//*************************************************
		if(object.getSource() == btnStop ){		btnStop.setIcon(new ImageIcon("img\\stop_hov.png"));	}
		if(object.getSource() == btnPlay ){		btnPlay.setIcon(new ImageIcon("img\\play_hov.png"));	}
		if(object.getSource() == btnPause ){	btnPause.setIcon(new ImageIcon("img\\pause_hov.png"));	}
		
		mouseDown=false;
		
	}

	@Override
	public void mouseDragged(MouseEvent obj) {
		//System.out.println("mouse dragged");
		mx = (int) Math.max(0, (obj.getX()/scaleOut));
		my = (int) Math.max(0, ((sy-obj.getY())/scaleOut));
		
		//System.out.println("mx "+mx);
		//System.out.println("my "+my);
		
		if(mDens && mx>2 && my>2 && ssx-2>mx && (ssy-2)>my){
			fs.dOld[mx][my]		=1;
			fs.dOld[mx-1][my]	=1;
			fs.dOld[mx+1][my]	=1;
			fs.dOld[mx][my-1]	=1;
			fs.dOld[mx][my+1]	=1;
			fs.dOld[mx-2][my]	=1;
			fs.dOld[mx+2][my]	=1;
			fs.dOld[mx][my-2]	=1;
			fs.dOld[mx][my+2]	=1;
			
			fs.dOld[mx-1][my+1]	=0.5f;
			fs.dOld[mx+1][my+1]	=0.5f;
			fs.dOld[mx-1][my-1]	=0.5f;
			fs.dOld[mx+1][my-1]	=0.5f;
		}
		
		if(mVel && mx>2 && my>2 && ssx-2>mx && (ssy-2)>my){
			
		fs.uOld[mx][my]		=mx-mxOld;
		fs.vOld[mx][my]		=my-myOld;
		
		fs.uOld[mx+1][my]	=mx-mxOld;
		fs.uOld[mx][my+1]	=mx-mxOld;
		fs.uOld[mx-1][my]	=mx-mxOld;
		fs.uOld[mx][my-1]	=mx-mxOld;
		
		fs.vOld[mx+1][my]	=my-myOld;
		fs.vOld[mx][my+1]	=my-myOld;
		fs.vOld[mx-1][my]	=my-myOld;
		fs.vOld[mx][my-1]	=my-myOld;
		
		myOld=my;
		mxOld=mx;
		}
		
	}

	@Override
	public void mouseMoved(MouseEvent obj) {
	
		
		mx = (int) Math.max(0,(obj.getX()/scaleOut));
		my = (int) Math.max(0,((sy-obj.getY())/scaleOut));
		
		if(dispDebug)refreshDebug(mx,my);
	}
	
	public float textToFloat(String value, float b0, float b1, float old, float scale){
		float x; 
		try {
		    x = Float.parseFloat(value);
		    
		} catch (NumberFormatException e) {
			x=old;
		    System.out.println("No valied number entered.");
		}
		
		if(x>b1)x=b1;
		if(x<b0)x=b0;
		x /= scale;
		return x;
	}
}
		