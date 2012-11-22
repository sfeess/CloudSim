
public class FluidSolver {
	
		
	//Gridsize (without +2 Boarder)
	int ssx, ssy;
	
	//Pixelsize
	int sx, sy;
	
	//density
	float[] d, dOld;
	
	//output
	float[] out;
	
	//velocity in pixels/seconds
	float[] u, uOld;
	float[] v, vOld;
	
	int size;
	
	// time in seconds
	float time;
	float dt;
	
	 /**Interpolation: 
     *  Search neighbor fields of given field position [xpos,ypos].
     *  Then interpolate value on xpos ypos
     **/
	public float interpolate(float xpos, float ypos, float[] f){
		
		xpos+=0.5;
		ypos+=0.5;
		
		// cellposition
		int x1 = (int) (xpos); 
		int x2 = (int) Math.abs(Math.floor(xpos)+(int)Math.signum(xpos%1-0.5));
		int y1 = (int) ypos;
		int y2 = (int) Math.abs(Math.floor(ypos)+(int)Math.signum(ypos%1-0.5));
		
		xpos-=0.5;
		ypos-=0.5;
		
		float dx = (float) Math.abs((xpos%1));
		float dy = (float) (1-Math.abs((ypos%1)));
		
		//Boarder Conditions
		if (x1<0.5){x1=0; x2=0;}
		if (x1>ssx+0.5){x1=ssx+1; x2=ssx+1;}
		if (y1<0.5){y1=0; y2=0; }
		if (y1>ssy+0.5){y1=ssy+1; y2=ssy+1;}
		
		
		//System.out.println(x1+" y"+y1);
		float value=		f[flin((int)x1,(int)y2)]*(1-dx)*(1-dy)+
							f[flin((int)x2,(int)y2)]*dx*(1-dy)+
							f[flin((int)x1,(int)y1)]*dy*(1-dx)+
							f[flin((int)x2,(int)y1)]*dx*dy;
		return value;
		
	}
	
	
	/**
	 * Evaluate: Scales field f to the size x, y;
	 * @param x = Targetsize X
	 * @param y = Targetsize Y
	 * @param f = SourceField to scale (dimensions=sx,sy) 
	 * @return interpolated field with dimensions x,y
	 */
	public float[] evaluate(float x, float y, float[] f){
		// field f hat die dimensions sx und sy  und soll auf x, y gemappt werden
		
		//mapped = result field 
		float[] mapped = new float[(int) (x*y)];
				
		//mapping in x richtung
		float xscale=((this.ssx+2)/x);
		float yscale=((this.ssy+2)/y);
		
		// gehe alle pixel durch
		for(int i=1; i<x-1; i++){
			for(int j=1; j<y-1; j++){
				//System.out.println("x"+i+" y"+j);				
				//transform Pixle positions in Cellspace 
				float xpos= xscale*i;
				float ypos= yscale*j;

				mapped[FluidViewer.plin(i,j)]= interpolate(xpos,ypos,f);	
			}
		}
		return mapped;
	}
	
	
	/** Step:
	 *  Execute Velocity and Density Solver
	 *  and update time and Output Field
	 **/
	public void step(){
		time=time+dt;
		solveVel();
		solveDens();
		
		out=d;
	}
	
	/** returns linearized ArrayIndex for Simulation Fields of position [i,j]
	 * 
	 */
	public int flin(int i, int j){		return ((i)+(this.ssx+2)*(j));}
	

	/**Setup the FluidSolver
	 * @param ssx = SimulationGrid size x (without boarder)
	 * @param ssy = SimulationGrid size y (without boarder)
	 * @param dt = Simulation time step in seconds
	 * @param sx = outputPixel size x
	 * @param sy = outputPixel size y
	 */
	public void setup(int ssx, int ssy, float dt, int sx, int sy){
		this.ssx = ssx;
		this.ssy = ssy;
		this.sx = sx;
		this.sy = sy;
		this.size = (ssx+2)*(ssy+2);
		this.dt = dt;
		reset();
	}
		
	/** Reset the FluidSolver => Set all fields to zero
	 * and set initial Field values
	 */
	public void reset(){
		u = uOld = v = vOld = d = dOld = out = new float[size];
		for(int x= 0; x<size; x++){
			u[x] = uOld[x] = v[x] = vOld[x] = d[x] = dOld[x] = out[x] = 0.0f;
		}
		//TEST=======================================
		//u = uOld = Field.circFieldU(ssx,ssy,1);
		v = vOld = Field.circFieldV(ssx,ssy,1);
		u = uOld = Field.circFieldU(ssx,ssy,1);
		//v = vOld = Field.constField(ssx,ssy,1F);

		d =dOld = Field.lineField(ssx,ssy);
		//System.out.println("d lenghhth="+d.length);
		out=d;
	}
	
	
	/**Swaps fields u and uOld */
	public void swapU(){float[] temp = u;	u=uOld;	uOld=temp;}
	/**Swaps fields v and vOld */
	public void swapV(){float[] temp = v;	v=vOld;	vOld=temp;}
	/**Swaps fields d and dOld */
	public void swapD(){float[] temp = d;	d=dOld;	dOld=temp;}
	
	
	
	
	
	public void addForce(){}
	
	/**
	 * Advects a field by the velocity field
	 * @param a = Field to advect
	 * @param aOld = Field to advect Old
	 * @param u = U-Velocity field
	 * @param v = V-Velocity field
	 */
	public void advect(float[] a, float[] aOld, float[] u, float[] v){
		
		float x,y,dtx,dty;
		
		//dtx=dt*ssx;
		//dty=dt*ssy;

		for (int i=1; i<=ssx; i++){
			for (int j=1; j<=ssy; j++){
				
				// gesuchte xKoordinate x = aktuell PunktX - dt * velocityU
				x = (float)i-dt*u[flin(i,j)];
				// gesuchte yKoordinate y = aktuell PunktY - dt * velocityV
				y = (float)j-dt*v[flin(i,j)];
								
				//Boarder Conditions
				if(x<0){x=0.5f;}
				if(y<0){y=0.5f;}
				if(x>ssx+1){y=ssx+0.5f;}
				if(y>ssy+1){y=ssy+0.5f;}
				
				float res = interpolate(x,y,aOld);
				a[flin(i,j)] = res;
			}
		}
		
	}
	
	
	public void diffuse(){}
	public void project(){}
	
	public void addSource(){}
	
	
	public void solveVel(){
		addForce();
		
		swapU();
		copy(u,uOld);
		//advect(u, u, v);
		
		swapV();
		copy(v,vOld);		
		//advect(v, u, v);
		
		diffuse();
		project();
		}
	
	public void solveDens(){
		
		addSource();
		diffuse();
		
		swapD();
		
		advect(d,dOld,u,v);
		//copy(d,dOld);
	}
	
	
public void copy(float[] _a, float[] _aOld){
		for (int i=0; i<ssx; i++){
			for (int j=0; j<ssy; j++){
				_a[flin(i,j)] = _aOld[flin(i,j)];
			}
		}
		
	}
	
	
	//
	
	

	
	
	 
	
}
