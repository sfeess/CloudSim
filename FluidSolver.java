
public class FluidSolver {

	//Gridsize (without +2 Boarder)
	int ssx, ssy;
	
	//Pixelsize
	int sx, sy;
	
	//Field length
	int size;
	
	//density
	float[] d, dOld;
	
	//output
	float[] out;
	
	//velocity in cells/sec
	float[] u, uOld;
	float[] v, vOld;
	
	// Vorticity
	float[] vorticity;
	
	// field to swap
	float[] temp;
	
	/** (nabla)q = pure Divergence component of the velocity field:  vel-(nabla)q = divergence free */
	float[] q;
	/** Divergence of one cell 	 */
	float[] div;
	
	// Diffusionrate & Viscosity 
	float diff;
	float visc;
	
	// time in seconds
	float time;
	float dt;

	int step;
	
	
	/**Setup the FluidSolver
	 * @param ssx = SimulationGrid size x (without boarder)
	 * @param ssy = SimulationGrid size y (without boarder)
	 * @param dt = Simulation time step in seconds
	 * @param sx = outputPixel size x
	 * @param sy = outputPixel size y
	 */
	public void setup(int ssx, int ssy, float dt, float scale){
		
		this.ssx = ssx;
		this.ssy = ssy;
		this.sx = (int) (ssx*scale);
		this.sy = (int) (ssy*scale);
		this.size = (ssx+2)*(ssy+2);
		this.dt = dt;
		reset();
	}
	
	
	/** Reset the FluidSolver => Set all fields to zero
	 * and set initial Field values
	 */
	public void reset(){
		
		step=0;
		
		u = new float[size];
		uOld = new float[size]; 
		v = new float[size];
		vOld = new float[size];
		d = new float[size];
		dOld = new float[size];
		out = new float[size];
		q = new float[size];
		div = new float[size];
		temp = new float[size];
		dOld = new float[size];
		vorticity = new float[size];
		
		for(int x= 0; x<size; x++){
			u[x] = uOld[x] = v[x] = vOld[x] = d[x] = dOld[x] = out[x] = q[x] = vorticity[x] = 0.0f;
		}
		
		//Field Init Setup=======================================
		v = Field.boxField(ssx,ssy,2f);
		u = Field.boxField(ssx,ssy,0f);
		
		
		d = Field.boxField(ssx,ssy,1f);
		
		diff = 0.0000f;
		visc = 0.000000001f;
		
		System.arraycopy(d,0,out,0,size);
	}
	
	
	/** Step:
	 *  Execute Velocity and Density Solver
	 *  and update time and Output Field
	 **/
	public void step(){
		time=time+dt;
		step++;
		solveVel();
		solveDens();
		
		System.arraycopy(d,0,out,0,size);
	}
	
	
	/**
	 * Velocity Solver
	 * addForce - viscDiffuse - Advect - Project 
	 */
	public void solveVel(){
		
		addSource(u,uOld);
		addSource(v,vOld);
		
		// add moving vel Source
		addSource(v,Field.smlBoxField(ssx,ssy,0.5f));
		addSource(u,Field.smlBoxField(ssx,ssy,(float) Math.sin(step*0.2)));
		
		vorticityConf(uOld, vOld, 0.2f);
		addSource(u,uOld);
		addSource(v,vOld);
		
		swapU(); swapV();
		copy(u,uOld); copy(v,vOld);	
		//diffuse(1, u, uOld, visc, dt);
		//diffuse(2, v, vOld, visc, dt);
		
		project(u,v,q,div);
		
		swapU(); swapV();
		//copy(u,uOld); copy(v,vOld);	
		advect(1, u, uOld, uOld, vOld);
		advect(2, v, vOld, uOld, vOld);
		
		project(u, v,q,div);
		for (int i = 0; i < size; i++) {uOld[i] = 0;vOld[i] = 0;}
		
		}
	
	
	/**
	 * Density Solver
	 * addDensity - Diffuse - Advect
	 */
	public void solveDens(){
		addSource(d, dOld);
		
		// add moving vel Source
		//addSource(d,Field.smlBoxField(ssx,ssy,0.5f));
		
		swapD();
		//copy(d,dOld); 
		diffuse(0, d, dOld, diff, dt);
		
		swapD();
		
		advect(0,d,dOld,u,v);
		//copy(d,dOld);
		
		for (int i = 0; i < size; i++) {dOld[i] = 0;}
    
	}
	
	
	/**
	 * Adds a field fOld to the field f 
	 * f += timestep * fOld 
	 * @param f0 = Field to add value to
	 * @param f1 = Field to be added
	 */
	public void addSource(float[] f0, float[] f1){
		for (int i = 0; i < size; i++){
	            f0[i] += dt * f1[i];}
	}
	
	
	/**
	 * Diffuse with the Gauss Seidel Relaxation
	 * x[i,j] = x0[i,j] + a * ( x[i-1,j] + x[i+1,j] + x[i,j-1] + x[i,j+1] ) / (1 + 4 * a)
	 */
	public void diffuse(int b, float[] f, float[] fOld, float diff, float dt ){
		
		int i, j, k; 
		float a=dt*diff*ssx*ssy; 
		
		//20 Relaxation iterations
		for ( k=0 ; k<20 ; k++ ) { 
			for ( i=1 ; i<=ssx ; i++ ) { 
				for ( j=1 ; j<=ssy ; j++ ) {    
					//Calc the average of cell + diff * neighbor cells
					f[flin(i,j)] = 	(fOld[flin(i,j)] + a*( 	f[flin(i-1,j)] + 
															f[flin(i+1,j)] + 
															f[flin(i,j-1)] + 
															f[flin(i,j+1)]  ) )
																				/(1+4*a); 
				} 
			} 
			setBounds(b,f); 
		} 
	}
	
	
	/**
	 * Advects a field by the velocity field
	 * @param b = BoundCond and Staggered Grid Type|| for: uVel=1(left celledge) ; vVel=2(top celledge) ; other(dens)(cellcenter)=0
	 * @param a = Field to advect
	 * @param aOld = Field to advect Old
	 * @param u = U-Velocity field
	 * @param v = V-Velocity field
	 */
	public void advect(int b, float[] a, float[] aOld, float[] u, float[] v){
		
		float x,y;

		// for density (or other central data)
		if(b==0){
			for (int i=1; i<=ssx; i++){
				for (int j=1; j<=ssy; j++){
					// x and y are in middle of the u/v valuepositions.
					x = (i+0.5f) - dt * ( u[flin(i,j)] + u[flin(i+1,j)] ) ;
					y = (j+0.5f) - dt * ( v[flin(i,j)] + v[flin(i,j+1)] );
					//Boarder Conditions
					if(x<0.0){		x=0.0f;}
					if(y<0.0){		y=0.0f;}
					if(x>ssx+0.5){	x=ssx+0.5f;}
					if(y>ssy+0.5){	y=ssy+0.5f;}
					// interpolate the value of old field
					a[flin(i,j)] = (float)interpolate(x-0.5f,y-0.5f,aOld);
				}
			}
		}
		
		// for u velocity (stored on left cell edge)
		if(b==1){
			for (int i=1; i<=ssx; i++){
				for (int j=1; j<=ssy; j++){
					//x is onValue - for y interpolate 4 surrounding v values
					x = i - dt * u[flin(i,j)];
					y = (j+0.5f) - dt * ( v[flin(i-1,j)] + v[flin(i,j)] + v[flin(i-1,j+1)] + v[flin(i,j+1)])/4;
					//Boarder Conditions
					if(x<0.5){		x=0.5f;}
					if(y<0.5){		y=0.5f;}
					if(x>ssx+0.5){	x=ssx+0.5f;}
					if(y>ssy+0.5){	y=ssy+0.5f;}
					// interpolate the value of old field
					a[flin(i,j)] = (float)interpolate(x,y-0.5f,aOld);
				}
			}
		}
		
		// for v velocity (stored on upper cell edge)
		if(b==2){
			for (int i=1; i<=ssx; i++){
				for (int j=1; j<=ssy; j++){
					// backtracked x&y coordinates xy = pointXY - dt * velocityUV
					// 1 -> position of v velocity
					// 2 -> interpolate velocity for Staggered grid 
					//__|==1==|__________|==============================2=====================================|
					x = (i+0.5f) - dt * ( u[flin(i,j-1)] + u[flin(i,j)] + u[flin(i+1,j-1)] + u[flin(i+1,j)])/4 ;
					y =  j -       dt *  v[flin(i,j)];
					//Boarder Conditions
					if(x<0.5){		x=0.5f;}
					if(y<0.5){		y=0.5f;}
					if(x>ssx+0.5){	x=ssx+0.5f;}
					if(y>ssy+0.5){	y=ssy+0.5f;}
					// interpolate the value of old field
					a[flin(i,j)] = (float)interpolate(x-0.5f,y,aOld);
				}
			}
		}
		setBounds(b,a);
	}


	/**
	 * Gauss Seidel Relaxation Solver
	 * Solves Equations in form  Ax=b with b=0 
	 * A  = Matrix
	 * a0 = matrix to store temp results
	 * 
	 */
	float[] linSolver(float[][] a, float[][] a0, float[] b){
		
		float u, l;
		
		// approximation x (at iteration step k) - xOld (at iteration step k-1)
		float[] x = new float[a.length];
		float[] xOld = new float[a.length];
		
		for(int i=0; i<x.length; i++ ){
			x[i]=0;  
			xOld[i]=0;  
		}
		
		for(int k=0; k<17; k++){
			// step through all lines 
			for(int i=0; i<a.length; i++){
				l=0;
				u=0;
				// step through all x in line
				for(int j=0; j<a.length; j++){
					// calc upper and lower triangle matrix
					if(j<i){	l+=a[i][j]*x[j];	}
					if(j>i){	u+=a[i][j]*xOld[j];	}
				}
				// calc new approximation
				x[i]= (1/a[i][i]) * (b[i]-l-u);
				
			}
			// stop when result is converged 
			if (x[0]==xOld[0]){	
				System.out.println("Ende nach "+k+" Iterationen");
				break;
				}
			// copy approximation x to xOld
			System.arraycopy(x,0, xOld, 0, a.length);
		}
		return x;
	}
	
	
	

	
	
	/**
	 * 
	 * 
	 * @param u = field to store projected u-vel field 
	 * @param v = field to store projected v-vel field 
	 * @param p = temporary pressure field  
	 * @param div = temporary divergence field  
	 */
	public void project(float[] u, float[] v,float[] q, float[] div){
		
		int i, j, k; 
		float h = 1.0f;  
		
		//calculate divergence field
		for ( i=1 ; i<=ssx ; i++ ) {   
			for ( j=1 ; j<=ssy ; j++ ) { 
				//old centered data
				//div[flin(i,j)] = 0.5f*h*(u[flin(i+1,j)]-u[flin(i-1,j)]+  v[flin(i,j+1)]-v[flin(i,j-1)]);    
				
				// div at cell center
				div[flin(i,j)] = 1f*h*(	u[flin(i+1,j)]	// u[flin(i+1,j)] 
										- 	u[flin(i,j)]	//u[flin(i-1,j)]
										+ 	v[flin(i,j+1)] //v[flin(i,j+1)]
										- 	v[flin(i,j)] 	//v[flin(i,j-1)]
										);  
				
				//initialize q for gauss seidel start value
				q[flin(i,j)] = 0; 
				
			}
		}
		setBounds(0, div);
		setBounds(0, q);
				
		//gauss seidel solve for q (i&j without boarders 0&ssxy+1)	
		for ( k=0 ; k<5 ; k++ ) {   
			for ( i=1 ; i<=ssx ; i++ ) {    
				for ( j=1 ; j<=ssy ; j++ ) { 
					q[flin(i,j)] = (-(div[flin(i,j)])
									+q[flin(i-1,j)]+q[flin(i+1,j)]
									+q[flin(i,j-1)]+q[flin(i,j+1)])/4;
					} 
			} 
		setBounds(0, q);  
		}
		
		// Subtract (nabla)q from u and v (makes vel divergence free)
		// u-(nabla)qx = u-(q right - q left)/2h
		// v-(nabla)qy = v-(q down - q up)/2h
		for ( i=1 ; i<=ssx ; i++ ) {   
			for ( j=1 ; j<=ssy ; j++ ) {  
				
				// central derivation? for u with i and i-1 ? or i-1 and i+1
				//u[flin(i,j)] -= (q[flin(i+1,j)]-q[flin(i-1,j)])/(2*h);    
				//v[flin(i,j)] -= (q[flin(i,j+1)]-q[flin(i,j-1)])/(2*h); 

				u[flin(i,j)] -= (q[flin(i+1,j)]-q[flin(i-1,j)])/(2*h);    
				v[flin(i,j)] -= (q[flin(i,j+1)]-q[flin(i,j-1)])/(2*h);   
			}  
		}
		setBounds(1, u);
		setBounds(2, v);
	} 
	
	/**
	 * Calculate vorticity at point i,j and store in vorticity field
	 * @param i
	 * @param j
	 * @return 
	 */
	public float vorticity(int i, int j){
		
		float du_dy =( 	(v[flin(i,j)] + v[flin(i+1,j)] + v[flin(i,j+1)] + v[flin(i+1,j+1)]) / 4 
				 	  - (v[flin(i,j)] + v[flin(i-1,j)] + v[flin(i,j+1)] + v[flin(i-1,j+1)]) / 4 ) / 2;
		
		float dv_dx = ( (u[flin(i,j+1)] + u[flin(i+1,j+1)] + u[flin(i,j)] + u[flin(i+1,j)]) / 4 
					  - (u[flin(i,j-1)] + u[flin(i+1,j-1)] + u[flin(i,j)] + u[flin(i+1,j)]) / 4 ) / 2;	
		
		
		//central diff of v in x direction - central diff of u in y direction
		return du_dy-dv_dx;
	}

	
	/** Vortex Confinement; increases the swirly motion that tends to disappear trough numerical dampening.
	 * 
	 * @param vF_u Velocity_u
	 * @param vF_v Velocity_v
	 * @param k	Confinement Multiplier
	 */
	public void vorticityConf(float[] vF_u, float[] vF_v, float k){

		float nab_nx, nab_ny, mag_n, nx, ny;
		

		//Calculate vorticity magnitude field = n
		for (int i=1; i<=ssx; i++){
			for (int j=1; j<=ssy; j++){
				vorticity[flin(i,j)]= Math.abs(vorticity(i,j));
			}
		}
		
		//Calculate vorticity Gradient N = (nabla n) / ||n||
		for (int i=1; i<=ssx; i++){
			for (int j=1; j<=ssy; j++){
				//calc nabla n for x and y by central difference
				nab_nx = (vorticity[flin(i+1,j)] - vorticity[flin(i-1,j)]) / 2;
				nab_ny = (vorticity[flin(i,j+1)] - vorticity[flin(i,j-1)]) / 2;
				
				//calculate magnitude of the vector (nab_nx , nab_ny)
				//add small value to prevent divide by zero
				mag_n = (float) Math.sqrt(nab_nx*nab_nx + nab_ny*nab_ny + 0.00000001f);
				
				// N = (nabla n) / ||n||
				nx=nab_nx/mag_n;
				ny=nab_ny/mag_n;
				
				// F = N x vorticity
				// F = Nx * vort_y - Ny * vort_x
				vF_u[flin(i,j)] = nx * vorticity(i,j)*k;
				vF_v[flin(i,j)] = - ny * vorticity(i,j)*k;
					
				
			}
		}
		
		
		
	}
		
	/** copy xOld to x	 */
	public void copy(float[] x, float[] xOld)	{System.arraycopy(xOld,0,x,0,xOld.length);}
	
	/**Swaps fields u and uOld */
	public void swapU(){temp = u; u=uOld;	uOld=temp;}
	/**Swaps fields v and vOld */
	public void swapV(){temp = v;	v=vOld;	vOld=temp;}
	/**Swaps fields d and dOld */
	public void swapD(){temp = d;	d=dOld;	dOld=temp;}
	
	/** returns linearized ArrayIndex for Simulation Fields of position [i,j] */
	public int flin(int i, int j){		return ((i)+(this.ssx+2)*(j));}
	
	
	/**
	 * Calc boundary values (bv) out of neighbor values (nv)
	 * @param b boundary type:	b=1:  bv=-nv (in x direction)
	 * 							b=2:  bv=-nv (in y direction)
	 * 							else: bv=nv 
	 * 							for velocity bounds b=1 or b=2 for pressure b=0
	 * @param f =field to correct bounds
	 */
	public void setBounds (int b, float[] f) { 
		// Border-column 0 and ssx+1 are being set
		for (int i=1 ; i<=ssy ; i++ ) {  
			f[flin(1, i)]    = b==1 ? -f[flin(2,i)]   : f[flin(1,i)];
			f[flin(0, i)]    = b==1 ? -f[flin(2,i)]   : f[flin(1,i)];   
			f[flin(ssx+1,i)] = b==1 ? -f[flin(ssx,i)] : f[flin(ssx,i)];   
		}
		
		// Border-line 0 and ssy+1 are being set   
		for (int i=1 ; i<=ssx ; i++ ) { 
			f[flin(i,1  )]   = b==2 ? -f[flin(i,2)]   : f[flin(i,1)]; 
			f[flin(i,0  )]   = b==2 ? -f[flin(i,2)]   : f[flin(i,1)]; 
			f[flin(i,ssy+1)] = b==2 ? -f[flin(i,ssy)] : f[flin(i,ssy)]; 
		} 
		//Corners - average between both direct neighbors
		f[flin(0  ,0  )] 	= 0.5f * (f[flin(1,0)]       + f[flin(0,1)]); 
		f[flin(0  ,ssy+1)] 	= 0.5f * (f[flin(1,ssy+1)]   + f[flin(0,ssy)]); 
		f[flin(ssx+1,0  )] 	= 0.5f * (f[flin(ssx,0)]     + f[flin(ssx+1,1)]); 
		f[flin(ssx+1,ssy+1)]= 0.5f * (f[flin(ssx,ssy+1)] + f[flin(ssx+1,ssy)]); 
		
	} 
	

	/**Interpolation: 
     *  Search neighbor fields of given field position [xpos,ypos].
     *  Then interpolate value on xpos ypos
     *  
     *  value not in center ==> values on integers (on gridlines)
     *  
     **/
	public float interpolate(float xpos, float ypos, float[] f){
		
		//Boarder Conditions
		if (xpos<0.5){xpos=0.5f;}
		if (xpos>ssx+0.5){xpos=ssx+0.5f;}
		if (ypos<0.5){ypos=0.5f;}
		if (ypos>ssy+0.5){ypos=ssy+0.5f;}
		
		// Sample positions
		int x1 = (int) xpos;
		int x2 = x1+1;
		int y1 = (int) ypos;
		int y2 = y1+1;
		
		// Distances
		float dx = xpos%1;
		float dy = ypos%1;
		
		// Interpolated Value
		float value = 	(1-dx) *( 	f[flin(x1, y1)]*(1-dy) 	+ 	f[flin(x1, y2)]*dy  	)
						+  dx  *(	f[flin(x2, y1)]*(1-dy)	+  	f[flin(x2, y2)]*dy  	);
		
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
		
		float[] mapped = new float[(int) (x*y)];
				
		float xscale=((this.ssx)/x);
		float yscale=((this.ssy)/y);
		
		for(int i=0; i<x; i++){
			for(int j=1; j<y; j++){
				// calc Sampleposition
				float xpos= xscale*i;
				float ypos= yscale*j;
				// interpolate at [xpos,ypos]
				mapped[FluidViewer.plin(i,j)]= interpolate(xpos,ypos,f);	
			}
		}
		return mapped;
	}
	
}
