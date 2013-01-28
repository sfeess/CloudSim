
public class FluidSolver {
	
	// Simulator quantities
	//************************************************************
	int sx, sy;					//Pixelsize
	int ssx, ssy;				// Gridsize (without +2 Boarder)
	int size;					// Field length
	int step;					// SimulationStep
	float time,dt;				// time and timestep in sec
	float scale;	
	
	// Fluid variables	
	//************************************************************
	float[][] u, uOld, v, vOld;	// velocity in cells/sec
	float[][] d, dOld;			// density
	float[][] vorticity;		// Vorticity
	float[][] temp;				// field to swap
	float diff; 				// Viscosity 
	float visc; 				// Diffusionrate
	float wind;					// horzontal wind
	
	// Cloud variables
	//************************************************************
	float[][] qc, qcOld;		// Condensed cloud water mixing ratio
	float[][] qv, qvOld;		// Water vapor mixing ratio
	float[][] qs;				// saturation vapor mixing ratio
	float[][] pt, ptOld;		// potential temperature
	float t0;					// ground temp
	float rd;					// gas constant in J/(kg K) 
	float maxAlt;				// upper boarder of altitude in simulation in meter
	float tlr;					// Temperature lapse rate in °C per 100meter
	float hum;					// humidity in percent 0-1
	float[] absT;				// absolute Temperature at altitude in K
	float [] absP;				// absolute Pressure at altitude in kPa
	
	
	// Output Fields
	//************************************************************
	int x1,x2,y1,y2;			// interpolation positions
	float dx,dy;				// interpolation distances
	float xpos, ypos;			// interpolation results
	
	
	
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
		this.scale = scale;
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
		
		
		// Fluid variables
		//************************************************************
		u = 	new float[ssx+2][ssy+2];
		uOld = 	new float[ssx+2][ssy+2];
		v = 	new float[ssx+2][ssy+2];
		vOld = 	new float[ssx+2][ssy+2];
		d = 	new float[ssx+2][ssy+2];
		dOld = 	new float[ssx+2][ssy+2];
		temp = 	new float[ssx+2][ssy+2];
		dOld = 	new float[ssx+2][ssy+2];
		vorticity = new float[ssx+2][ssy+2];
		
		
		// Cloud variables
		//************************************************************
		qc = 	new float[ssx+2][ssy+2];
		qcOld = new float[ssx+2][ssy+2];
		qv = 	new float[ssx+2][ssy+2];
		qvOld = new float[ssx+2][ssy+2];
		qs = 	new float[ssx+2][ssy+2];
		pt = 	new float[ssx+2][ssy+2];
		ptOld = new float[ssx+2][ssy+2];
		absT = 	new float[ssy+2];
		absP = 	new float[ssy+2];
		
		
		
		// Set Fields to Zero
		//************************************************************
		for(int x= 0; x<ssx+2; x++){
			for(int y= 0; y<ssy+2; y++){
				u[x][y] = uOld[x][y] = v[x][y] = vOld[x][y] = d[x][y] = dOld[x][y] = vorticity[x][y] = qc[x][y] = 0.0f;  // = out[x]
			}
		}
		// Initialize Fluid 
		//************************************************************
		wind = 0.0f;
		diff = 0.0000f;
		visc = 0.000000001f;
		
		//v = Field.boxField(ssx,ssy,0f);
		//u = Field.constField(ssx,ssy,0f);
		//d = Field.boxField(ssx,ssy,1f);
		
		
		// Initialize Cloud constants
		//************************************************************
		rd = 287; 		// specific gas constant

		// user defined
		maxAlt = 5000;	// altitude in meter on top of sim grid
		tlr = 0.6f; 	// Kelvin per 100 meter between 0.55 and 0.99
		tlr /= 100; 	// Kelvin per 1 meter
		t0 = 295;		// temp on ground
		
		hum = 0.5f;		// humidty
		
		// Initialize absolute Temp lookup
		//************************************************************
		for(int y= 0; y<ssy+2; y++){
			// Ground Temp - (altitude / 100)* tempLapseRate
			float alt= ((float)y/(float)ssy)*maxAlt;
			absT[y] = t0- alt*tlr;
		}
		
		// Initialize absolute Pressure lookup in kPa
		//************************************************************
		for(int y= 0; y<ssy+2; y++){
			float alt= ((float)y/(float)ssy)*maxAlt;
			absP[y] = (float) (10* Math.pow(( 1- (alt*tlr/t0) ),9.81/(tlr*rd))) ;
		}
		
		// Initialize pot temp
		//************************************************************
		for(int i= 0; i<ssx+2; i++){
			for(int j= 0; j<ssy+2; j++){
				pt[i][j] = (float) (absT[j]*Math.pow((100/absP[j]),0.286));  
				
			}
		}
		
		// Initialize Saturation vapor mixing ratio   and water vapor mixing ratio
		//************************************************************
		//T in celsius
		for(int i= 0; i<ssx+2; i++){
			for(int j= 0; j<ssy+2; j++){
					// temp in °C and p in Pa
					qs[i][j] = 		(float) Math.pow(380/(absP[j]*1000),(17.67*(absT[j]-273.15))/(absT[j]-273.15+243.5));
					qv[i][j] = 		qs[i][j] * hum;
					qvOld[i][j] = 	qs[i][j] * hum;
			}
		}
		
		
		
		//System.arraycopy(d,0,out,0,size);
	}
	
	
	/** Step:
	 *  Execute Velocity and Density Solver
	 *  and update time and Output Field
	 **/
	public void step(){
		time=time+dt;
		step++;
		solveVel();
		solveClouds();
		solveDens();
		//System.arraycopy(d,0,out,0,size);
	}
	
	
	/**
	 * Velocity Solver
	 * addForce - viscDiffuse - Advect - Project 
	 */
	public void solveVel(){
		
		addSource(u,uOld);
		addSource(v,vOld);
		
		//addSource(v,Field.noiseEmitField(ssx,ssy,0.51f,4f,time*0.1f,1f));
		// add moving vel Source
		//addSource(v,Field.smlBoxField(ssx,ssy,0.5f));
		//addSource(u,Field.smlBoxField(ssx,ssy,(float) Math.sin(step*0.2)));
		
		
		
		//vorticityConf(uOld, vOld, 0.1f);
		//addSource(u,uOld);
		//addSource(v,vOld);
		
		buoyancy(vOld, -0.000001f);
		addSource(v,vOld);
		
		swapU(); swapV();
		copy(u,uOld); copy(v,vOld);	
		//diffuse(1, u, uOld, visc, dt);
		//diffuse(2, v, vOld, visc, dt);
		project(u,v,uOld,vOld);
		
		swapU(); swapV();
		//copy(u,uOld); copy(v,vOld);	
		advect(1, u, uOld, uOld, vOld);
		advect(2, v, vOld, uOld, vOld);
		
		
		//addSource(qv,Field.noiseEmitField(ssx,ssy,0.51f,10f,time,(float)Math.abs(Math.sin(time))*0.01f));
		//Advect vapor, condensed water and pot_Temp
		swapQv(); swapQc(); swapPt();
		advect(3, qv, qvOld, u, v);
		advect(0, qc, qcOld, u, v);
		advect(0, pt, ptOld, u, v);
		
		waterCont();
		
		//qv= Field.noiseEmitField(ssx,ssy,0.51f,4f,time*1.1f,1f);
		
		project(u,v,uOld,vOld);
		
		// reset uOld vOld
		for (int i = 0; i < ssx+2; i++) {
			for (int j = 0; j < ssy+2; j++){
				uOld[i][j] = 0;vOld[i][j] = 0;
			}
		}
		
		
		}
	
	public void solveClouds(){
	
	}
	
	/**
	 * Density Solver
	 * addDensity - Diffuse - Advect
	 */
	public void solveDens(){
		addSource(d, dOld);
		
		// add moving vel Source
		//addSource(d,Field.noiseEmitField(ssx,ssy,0.51f,10f,time,1f));
		
		swapD();
		copy(d,dOld); 
		//diffuse(0, d, dOld, diff, dt);
		
		swapD();
		
		advect(0,d,dOld,u,v);
		//copy(d,dOld);
		
		
		for (int i = 0; i < ssx+2; i++) {
			for (int j = 0; j < ssy+2; j++){
				dOld[i][j] = 0;
			}
		}
		
	}
	
	
	/**
	 * Adds a field fOld to the field f 
	 * f += timestep * fOld 
	 * @param f0 = Field to add value to
	 * @param f1 = Field to be added
	 */
	public void addSource(float[][] f0, float[][] f1){
		
		for (int i = 0; i < ssx+2; i++) {
			for (int j = 0; j < ssy+2; j++){
				f0[i][j] += dt * f1[i][j];
			}
		}
		
	}
	
	
	/**
	 * Diffuse with the Gauss Seidel Relaxation
	 * x[i,j] = x0[i,j] + a * ( x[i-1,j] + x[i+1,j] + x[i,j-1] + x[i,j+1] ) / (1 + 4 * a)
	 */
	public void diffuse(int b, float[][] f, float[][] fOld, float diff, float dt ){
		
		int i, j, k; 
		float a=dt*diff*ssx*ssy; 
		
		//20 Relaxation iterations
		for ( k=0 ; k<20 ; k++ ) { 
			for ( i=1 ; i<=ssx ; i++ ) { 
				for ( j=1 ; j<=ssy ; j++ ) {    
					//Calc the average of cell + diff * neighbor cells
					f[i][j] = 	(fOld[i][j] + a*( 			f[i-1][j] + 
															f[i+1][j] + 
															f[i][j-1] + 
															f[i][j+1]  ) )
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
	public void advect(int b, float[][] a, float[][] aOld, float[][] u, float[][] v){
		
		float x,y;

		// for density (or other central data)
		if(b==0){
			for (int i=1; i<=ssx; i++){
				for (int j=1; j<=ssy; j++){
					// x and y are in middle of the u/v valuepositions.
					x = (i+0.5f) - dt * ( u[i][j] + u[i+1][j] ) ;
					y = (j+0.5f) - dt * ( v[i][j] + v[i][j+1] );
					//Boarder Conditions
					if(x<0.0){		x=0.0f;}
					if(y<0.0){		y=0.0f;}
					if(x>ssx+0.5){	x=ssx+0.5f;}
					if(y>ssy+0.5){	y=ssy+0.5f;}
					// interpolate the value of old field
					a[i][j] = (float)interpolate(x-0.5f,y-0.5f,aOld);
				}
			}
			setBounds(b,a);
		}
		
		// for u velocity (stored on left cell edge)
		if(b==1){
			for (int i=1; i<=ssx; i++){
				for (int j=1; j<=ssy; j++){
					//x is onValue - for y interpolate 4 surrounding v values
					x = i - dt * u[i][j];
					y = (j+0.5f) - dt * ( v[i-1][j] + v[i][j] + v[i-1][j+1] + v[i][j+1])/4;
					//Boarder Conditions
					if(x<0.5){		x=0.5f;}
					if(y<0.5){		y=0.5f;}
					if(x>ssx+0.5){	x=ssx+0.5f;}
					if(y>ssy+0.5){	y=ssy+0.5f;}
					// interpolate the value of old field
					a[i][j] = (float)interpolate(x,y-0.5f,aOld);
				}
			}
			setBounds(b,a);
		}
		
		// for v velocity (stored on upper cell edge)
		if(b==2){
			for (int i=1; i<=ssx; i++){
				for (int j=1; j<=ssy; j++){
					// backtracked x&y coordinates xy = pointXY - dt * velocityUV
					// 1 -> position of v velocity
					// 2 -> interpolate velocity for Staggered grid 
					//__|==1==|__________|==============================2=====================================|
					x = (i+0.5f) - dt * ( u[i][j-1] + u[i][j] + u[i+1][j-1] + u[i+1][j])/4 ;
					y =  j -       dt *  v[i][j];
					//Boarder Conditions
					if(x<0.5){		x=0.5f;}
					if(y<0.5){		y=0.5f;}
					if(x>ssx+0.5){	x=ssx+0.5f;}
					if(y>ssy+0.5){	y=ssy+0.5f;}
					// interpolate the value of old field
					a[i][j] = (float)interpolate(x-0.5f,y,aOld);
				}
			}
			setBounds(b,a);
		}
		
		// for qv periodic sides
		if(b==3){
			for (int i=0; i<=ssx; i++){
				for (int j=0; j<=ssy; j++){
					// x and y are in middle of the u/v valuepositions.
					x = (i+0.5f) - dt * ( u[i][j] + u[i+1][j] ) ;
					y = (j+0.5f) - dt * ( v[i][j] + v[i][j+1] );
					//System.out.println(i+"     "+x);
					//Boarder Conditions
					if(x<0.5){	
						//System.out.println("smaller 1");
						x=ssx+x;
					}
					if(y<0.0){		y=0.0f;}
					if(x>ssx+1.5){	
						System.out.println("greater ssx+1");
						x=ssx-x;}
					if(y>ssy+0.5){	y=ssy+0.5f;}
					// interpolate the value of old field
					a[i][j] = (float)interpolate(x-0.5f,y-0.5f,aOld);
				}
			}
			//setBounds(b,a);
		}
		
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
	public void project(float[][] u, float[][] v,float[][] q, float[][] div){
		
		int i, j, k; 
		float h = 1.0f;  
		
		//calculate divergence field
		for ( i=1 ; i<=ssx ; i++ ) {   
			for ( j=1 ; j<=ssy ; j++ ) { 
				//old centered data
				//div[flin(i,j)] = 0.5f*h*(u[flin(i+1,j)]-u[flin(i-1,j)]+  v[flin(i,j+1)]-v[flin(i,j-1)]);    
				
				// div at cell center
				div[i][j] = 1f*h*(			u[i+1][j]	// u[flin(i+1,j)] 
										- 	u[i][j]		//u[flin(i-1,j)]
										+ 	v[i][j+1] 	//v[flin(i,j+1)]
										- 	v[i][j]		//v[flin(i,j-1)]
										);  
				
				//initialize q for gauss seidel start value
				q[i][j] = 0; 
				
			}
		}
		setBounds(0, div);
		setBounds(0, q);
				
		//gauss seidel solve for q (i&j without boarders 0&ssxy+1)	
		for ( k=0 ; k<5 ; k++ ) {   
			for ( i=1 ; i<=ssx ; i++ ) {    
				for ( j=1 ; j<=ssy ; j++ ) { 
					q[i][j] = 	(-(div[i][j])
									+q[i-1][j]+q[i+1][j]
									+q[i][j-1]+q[i][j+1])/4;
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

				u[i][j] -= (q[i+1][j]	-	q[i-1][j])/(2*h);    
				v[i][j] -= (q[i][j+1]	-	q[i][j-1])/(2*h);   
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
	
		float dv_dx =(	(v[i][j] + v[i+1][j] + v[i][j+1] + v[i+1][j+1]) / 4 
			 	  - (v[i][j] + v[i-1][j] + v[i][j+1] + v[i-1][j+1]) / 4 ) * 0.5f;
	
		float du_dy = ( (u[i][j+1] + u[i+1][j+1] + u[i][j] + u[i+1][j]) / 4 
				  - (u[i][j-1] + u[i+1][j-1] + u[i][j] + u[i+1][j]) / 4 ) *0.5f;	
		
		
		//central diff of v in x direction - central diff of u in y direction
		return du_dy - dv_dx;
	}

	
	/** Vortex Confinement; increases the swirly motion that tends to disappear trough numerical dampening.
	 * 
	 * @param vF_u Velocity_u
	 * @param vF_v Velocity_v
	 * @param k	Confinement Multiplier
	 */
	public void vorticityConf(float[][] vF_u, float[][] vF_v, float k){

		float nab_nx, nab_ny, mag_n, nx, ny;
		

		//Calculate vorticity magnitude field = n
		for (int i=1; i<=ssx; i++){
			for (int j=1; j<=ssy; j++){
				vorticity[i][j]= Math.abs(vorticity(i,j));
			}
		}
		
		//Calculate vorticity Gradient N = (nabla n) / ||n||
		for (int i=2; i<ssx; i++){
			for (int j=2; j<ssy; j++){
				//calc nabla n for x and y by central difference
				nab_nx = (vorticity[i+1][j] - vorticity[i-1][j]) / 2;
				nab_ny = (vorticity[i][j+1] - vorticity[i][j-1]) / 2;
				
				//calculate magnitude of the vector (nab_nx , nab_ny)
				//add small value to prevent divide by zero
				mag_n = (float) Math.sqrt(nab_nx*nab_nx + nab_ny*nab_ny) + 0.00000001f;
				
				// N = (nabla n) / ||n||
				nx=nab_nx/mag_n;
				ny=nab_ny/mag_n;
				
				// F = N x vorticity
				// F = Nx * vort_y - Ny * vort_x
				vF_u[i][j] = - ny * (vorticity(i,j)+vorticity(i-1,j))*0.5f *k;
				vF_v[i][j] = nx * (vorticity(i,j)+vorticity(i,j-1))*0.5f*k;
					
				
			}
		}
		
		
		
	}
	
	public void initClouds(){
		
		
	}
	
	// Update the potential temperature according to condesation
	// Due to condensation latent energy is released in form of heat. -> change in pot temp
	public void thermoDyn(){
		
		
		
		
	}
	
	public void buoyancy(float[][] f, float k){
		float vpt;
		
		for (int i=1; i<=ssx; i++){
			for (int j=1; j<=ssy; j++){
				vpt = pt[i][j] * ( 1 + 0.61f * qv[i][j]);
				f[i][j] = k*( ( (vpt-295) / 295 ) - 9.81f * qc[i][j] );
			}
		}
	}
	
	
	public void waterCont(){
		
		float d_qv,T,qs;
		
		
		for (int i=1; i<=ssx; i++){
			for (int j=1; j<=ssy; j++){
				//alt=((float)j/(float)ssy)*maxAlt;
				
				//compute 	p(alt) in kPa
				// 			p(alt) = p0*(1- (zL/T0) )^(g/L Rd) 
				//			with alt= altitude in m, g = 9.81(m/s^2)
				//			p0 and T0 = p and T at the base altitude
				//			typically: p0 = 10kP and T0=280-310K
				// 			L = Lapse rate in °K or °C per meter
				//			Rd = ideal gas constant ~ 287 J/(kg K)
				
				//p = (float) Math.pow(10*(1-(alt*tlr/t0)),(9.81/(tlr*rd)));
				//
				//compute 	T = pt[i,j]/( (^p/p)^k  )     
				// 			with  ^p = 100kPa   k = ~0.286
				//			T = pt[i,j]/( (100/p)^0.286)
				T = (float) (pt[i][j]/Math.pow((100/absP[j]),0.286));
				// conversion from Kelvin to Celsius
				T = T-273.15f;
				//
				//compute 	qs= (380.16/p)exp((17.67*T)/(T+243.5))
				// with T in °C and P in Pa
				//
				// qs nicht als feld
				qs = (float) Math.pow((380.16/(absP[j]*1000)),((17.67*T)/(T+243.5)));
				
				d_qv  = Math.min(qs - qv[i][j],qc[i][j]);
				
				qv[i][j] = qv[i][j] + d_qv;
				qc[i][j] = qc[i][j] - d_qv;
				
				
				//update potential Temperature ( Thermodynamics Equation )
				//delta p = -L/(cp*PI)*(delta C)
				// L  = constant for latent heat released    	2.501 J/kg
				// cp = specific heat capacity of dry air 		1005 J/(kg K)
				// C  = condensation rate = condensation per evaporation
				//	  = - min(qvs-qv, qc)  siehe waterCont = d_qv
				// PI = Exner Function = T/pt
				//_______________-L_____/___cp___*________PI___________*____________C_________________________________
				T += 273.15f;
				pt[i][j] += -2.501 / ( 1005 * (T/pt[i][j]) ) * (-d_qv);
				
			}
		}
		
	}
	
	
	
	
	/** copy xOld to x	 */
	//
	public void copy(float[][] x, float[][] xOld)	
	{
		for (int i = 0; i < ssx+2; i++) {System.arraycopy(xOld[i], 0, x[i], 0, xOld[0].length);}
	}
	
	/**Swaps fields u and uOld */
	public void swapU(){temp = u; 	u=uOld;		uOld=temp;}
	/**Swaps fields v and vOld */
	public void swapV(){temp = v;	v=vOld;		vOld=temp;}
	/**Swaps fields d and dOld */
	public void swapD(){temp = d;	d=dOld;		dOld=temp;}
	/**Swaps fields pt and ptOld */
	public void swapPt(){temp = pt;	pt=ptOld;	ptOld=temp;}
	/**Swaps fields qc and qcOld */
	public void swapQc(){temp = qc;	qc=qcOld;	qcOld=temp;}
	/**Swaps fields qv and qvOld */
	public void swapQv(){temp = qv;	qv=qvOld;	qvOld=temp;}
	
	
	/** returns linearized ArrayIndex for Simulation Fields of position [i,j] */
	//public int flin(int i, int j){		return ((i)+(this.ssx+2)*(j));}
	
	
	/**
	 * Calc boundary values (bv) out of neighbor values (nv)
	 * @param b boundary type:	b=1:  bv=-nv (in x direction)
	 * 							b=2:  bv=-nv (in y direction)
	 * 							b=3:  periodic sidewalls
	 * 							else: bv=nv 
	 * 							for velocity bounds b=1 or b=2 for pressure b=0
	 * @param f =field to correct bounds
	 */
	public void setBoundsOriginal (int b, float[][] f) { 
		
		
		// Border-column 0 and ssx+1 are being set
		for (int i=1 ; i<=ssy ; i++ ) {  
			f[1][i]   	= b==1 ? -f[2][i]    : 	f[1][i] ;
			f[0][i]    	= b==1 ? -f[2][i]    : 	f[1][i] ;   
			f[ssx+1][i] = b==1 ? -f[ssx][i]  : 	f[ssx][i] ;   
		}
		
		// Border-line 0 and ssy+1 are being set   
		for (int i=1 ; i<=ssx ; i++ ) { 
			f[i][1]  	= b==2 ? -f[i][2]    : 	f[i][1]; 
			f[i][0]  	= b==2 ? -f[i][2]    : 	f[i][1] ; 
			f[i][ssy+1] = b==2 ? -f[i][ssy]  : 	f[i][ssy] ; 
		} 
		//Corners - average between both direct neighbors
		f[0][0]	= 0.5f * (f[1][0]      + f[0][1]); 
		f[0][ssy+1] 	= 0.5f * (f[1][ssy+1]   + f[0][ssy]); 
		f[ssx+1][0] 	= 0.5f * (f[ssx][0]     + f[ssx+1][1]); 
		f[ssx+1][ssy+1]	= 0.5f * (f[ssx][ssy+1] + f[ssx+1][ssy]); 
		
	} 
	
	
	public void setBounds (int b, float[][] f){
		
		// b=0 central data neuman boundary
		if (b==0){
			for (int i=1 ; i<=ssx ; i++ ) { 
				f[i][1]   = f[i][1]; 
				f[i][0]   = f[i][1]; 
				f[i][ssy+1] = f[i][ssy]; 
				//f[flin(i,0  )]   = 0.3f; 
				//f[flin(i,ssy+1)] = 0.3f; 
			}
			for (int i=1 ; i<=ssy ; i++ ) { 
				f[0][i]    = f[1][i];   
				f[ssx+1][i] = f[ssx][i];
				
			}
		}
		
		// b=1 u velocity
		// bottom 	= noslip
		// top 		= free slip 
		// sides 	= user defined wind
		if (b==1){
			for (int i=1 ; i<=ssy ; i++ ) { 
				f[1][i]    = -f[2][i];
				f[0][i]    = wind;   
				f[ssx+1][i] = wind; 
			}
			for (int i=1 ; i<=ssx ; i++ ) { 
				f[i][1]   =  f[i][1]; 
				f[i][0]   =  f[i][1]; 
				f[i][ssy+1] =  f[i][ssy]; 
			}
		}
		
		// b=2 v
		// bottom 	= noslip
		// top 		= free slip 
		// sides 	= zero
		if (b==2){
			for (int i=1 ; i<=ssy ; i++ ) { 
				f[0][i]     = 0;   
				f[ssx+1][i] = 0;  
			}
			for (int i=1 ; i<=ssx ; i++ ) { 
				f[i][1]   = -f[i][2]; 
				f[i][0]   = -f[i][2]; 
				f[i][ssy+1] = -f[i][ssy]; 
			}
		}
		
		// b=3 central data no slip
		if (b==3){
			
		}
		
	}

	/**Interpolation: 
     *  Search neighbor fields of given field position [xpos,ypos].
     *  Then interpolate value on xpos ypos
     *  
     *  value not in center ==> values on integers (on gridlines)
     *  
     **/
	public float interpolate(float xpos, float ypos, float[][] f){
		
		//Boarder Conditions
		if (xpos<0.5){xpos=0.5f;}
		if (xpos>ssx+0.5){xpos=ssx+0.5f;}
		if (ypos<0.5){ypos=0.5f;}
		if (ypos>ssy+0.5){ypos=ssy+0.5f;}
		
		// Sample positions
		x1 = (int) xpos;
		x2 = x1+1;
		y1 = (int) ypos;
		y2 = y1+1;
		
		// Distances
		dx = xpos%1;
		dy = ypos%1;
		
		// Interpolated Value
		return  	(1-dx) *( 	f[x1][y1]*(1-dy) 	+ 	f[x1][y2]*dy  	)
					+  dx  *(	f[x2][y1]*(1-dy)	+  	f[x2][y2]*dy  	);
		

		
	}
	
	
	/**
	 * Evaluate: Scales field f to the size x, y;
	 * @param x = Targetsize X
	 * @param y = Targetsize Y
	 * @param f = SourceField to scale (dimensions=sx,sy) 
	 * @return interpolated field with dimensions x,y
	 */
	public float[][] evaluate(float x, float y, float[][] f){
		
		float[][] mapped = new float[(int) x][(int) y];
				
		//float xscale=((this.ssx)/x);
		//float yscale=((this.ssy)/y);
		
		for(int i=1; i<x; i++){
			for(int j=1; j<y; j++){
				// calc Sampleposition
				//float xpos= xscale*i;
				//float ypos= yscale*j;
				// interpolate at [xpos,ypos]
				mapped[i][j]= interpolate(1/scale*i,1/scale*j,f);	
			}
		}
		return mapped;
	}
	
}
