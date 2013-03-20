

public class FluidSolver {
	
	float qc_max=0;
	int rkSteps=1;
	 
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
	float wind;					// horizontal wind
		
	
	// Cloud variables
	//************************************************************
	float[][] qc, qcOld;		// Condensed cloud water mixing ratio
	float[][] qv, qvOld;		// Water vapor mixing ratio
	float qs;					// saturation vapor mixing ratio
	float[][] pt, ptOld;		// potential temperature
	// Constants
	float tlr;					// Temperature lapse rate in °C per 100meter
	float hum;					// humidity in percent 0-1
	float maxAlt;				// upper boarder of altitude in simulation in meter
	float t0;					// ground temp
	float vort, buoyancy,rd,p0,grav,lh,cp,exner;
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
		init();
	}
	
	
	/** Reset the FluidSolver => Set all fields to zero
	 * and set initial Field values
	 */
	public void init(){
		
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
		pt = 	new float[ssx+2][ssy+2];
		ptOld = new float[ssx+2][ssy+2];
		absT = 	new float[ssy+2];
		absP = 	new float[ssy+2];
		
		
		// Initialize Fluid 
		//************************************************************
		wind = 0.0f;
		diff = 0.0000f;
		visc = 0.000000000f;
		
		
		
		// Initialize Cloud constants
		//************************************************************
		rd = 		287; 		// specific gas constant for dry air
		p0 = 		100;		// pressure at sea level (kPa)
		grav = 		9.81f;		// gravitational acceleration (m/s²)
		lh = 		2501;		// Latent heat of vaporization of water (J/kg)
		cp = 		1005;		// specific heat cpapcity J/(kg K)
		
		
		// user defined
		maxAlt = 	3000;		// altitude in meter on top of sim grid
		tlr = 		0.9f; 		// Kelvin per 100 meter between 0.55 and 0.99
			tlr /= 	100; 		// Kelvin per 1 meter
		t0 = 		295;		// temp on ground in Kelvin
		hum = 		0.6f;		// humidty
		buoyancy =  0.8f;
		vort = 		0.2f;
		
		reset();
	}	
		
		
	public void reset(){		
		
		step=0;
		
		// Set Fields to Zero
		//************************************************************
		for(int x= 0; x<ssx+2; x++){
			for(int y= 0; y<ssy+2; y++){
				u[x][y] = uOld[x][y] = v[x][y] = vOld[x][y] = d[x][y] = dOld[x][y] = vorticity[x][y] = qc[x][y] = 0.0f;  // = out[x]
			}
		}
		
		//v = Field.circFieldV(ssx,ssy,1f);
		//u = Field.constField(ssx,ssy,1f);
		//d = Field.boxField(ssx,ssy,1f);
		
		// Initialize absolute Temp lookup
		//************************************************************
		for(int y= 0; y<ssy+2; y++){
			// Ground Temp - (altitude / 100)* tempLapseRate
			float alt = ( (float) y / (float) ssy ) * maxAlt;
			absT[y] = t0 - alt * tlr;
		}
		
		// Initialize absolute Pressure lookup in kPa
		//************************************************************
		for(int y= 0; y<ssy+2; y++){
			float alt = ( (float)y / (float)ssy ) * maxAlt;
			absP[y] = (float) (p0* Math.pow(( 1- ( (alt*tlr)/t0 ) ),(grav/(tlr*rd)) )) ;
		}
		
		// Initialize pot temp
		//************************************************************
		for(int i= 0; i<ssx+2; i++){
			for(int j= 0; j<ssy+2; j++){
				//					K                   kPa/kPa
				pt[i][j] = (float) (absT[j] * ( Math.pow( (p0/absP[j]) , 0.286)));  
				ptOld[i][j] = pt[i][j]; //absT[j] * (float) ( Math.pow( (p0 / absP[j]) , 0.286)); 
			}
		}
		
		// Initialize Saturation vapor mixing ratio   and water vapor mixing ratio
		//************************************************************
		//T in celsius
		for(int i= 0; i<ssx+2; i++){
			for(int j= 0; j<ssy+2; j++){
					// temp in °C and p in Pa
					qs = 			(float) (  (380/(absP[j]*1000)  ) * Math.exp( (17.67*(absT[j]-273.15)) / (absT[j]-273.15+243.5))) ;
					qv[i][j] = 		qs * hum;
					qvOld[i][j] = 	qs * hum;
			}
		}
		
		
		
		
	}
	
	
	/** Step:
	 *  Execute Velocity and Density Solver
	 *  and update time and Output Field
	 **/
	public void step(){
		time=time+dt;
		step++;
		solveVel();
		//solveDens();
	}
	
	public void debugLine(String a){
		for(int j= 31; j>30&&j<60; j++){ 
			//System.out.println(a+" point"+j+"    value:"+d[20][j]);
		}
	}
	
	/**
	 * Velocity Solver
	 * addForce - viscDiffuse - Advect - Project 
	 */
	public void solveVel(){
		
		//addSource(qv,qvOld);
		
		//addSource(u,uOld);
		//addSource(v,vOld);
		
		//add moving vel Source
		//addSource(v,Field.smlBoxField(ssx,ssy,0.5f));
		//addSource(v,Field.smlBoxField(ssx,ssy,-Math.min(0,(float) Math.sin(step*0.1) ) ) );
		
		vorticityConf(uOld, vOld, vort);
		addSource(u,uOld);
		addSource(v,vOld);
		
		project(u,v,uOld,vOld);
		
		buoyancy(vOld, buoyancy);
		addSource(v,vOld);
		
		swapQv(); swapQc(); swapPt();
		//copy(qv,qvOld); copy(qc,qcOld);copy(pt,ptOld);
		advect(4, qv, qvOld, u, v);
		advect(5, qc, qcOld, u, v);
		advect(3, pt, ptOld, u, v);
		
		swapU(); swapV();
		//copy(u,uOld); copy(v,vOld);	
		advect(1, u, uOld, uOld, vOld);
		advect(2, v, vOld, uOld, vOld);
		
		waterCont();
		
		//diffuse(1, u, uOld, visc, dt);
		//diffuse(2, v, vOld, visc, dt);
		
		project(u,v,uOld,vOld);
		
		
		
		//addSource(qv,Field.noiseEmitField(ssx,ssy,0.51f,10f,time,(float)Math.abs(Math.sin(time))*0.01f));
		
		// reset uOld vOld
		for (int i = 0; i < ssx+2; i++) {
			for (int j = 0; j < ssy+2; j++){
				uOld[i][j] = 0;vOld[i][j] = 0; qvOld[i][j] = 0;
			}
		}
		
		
		}
	
	
	
	/**
	 * Density Solver
	 * addDensity - Diffuse - Advect
	 */
	
	
	public void solveDens(){
		
		//addSource(d, dOld);
		//debugLine("after sourceAdd");
		// add moving vel Source
		//addSource(d,Field.noiseEmitField(ssx,ssy,0.51f,10f,time,1f));
		//addSource(d,Field.smlBoxField(ssx,ssy,0.0415f));
		
		//swapD();
		//debugLine("after swap");
		//copy(d,dOld); 
		//debugLine("after copy");
		//diffuse(0, d, dOld, diff, dt);
		
		swapD();
		//debugLine("after swap");
		
		//advect(0,d,dOld,u,v);
		advect(5,d,dOld,u,v);
		//debugLine("after advection");
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
	 * @param b =  0=central Neumann || 1=uVel || 2=vVel || 3=potTemp || 4=waterVapor || 5=cloudWater
	 * @param a = Field to advect
	 * @param aOld = Field to advect Old
	 * @param u = U-Velocity field
	 * @param v = V-Velocity field
	 */
	public void advect(int b, float[][] a, float[][] aOld, float[][] u, float[][] v){
		
		float x, y, offset_x, offset_y;
		
		// offset to left edge
		if(b==1){
			offset_x = 0;
			offset_y = 0.5f;
		}
		// offset to top edge
		else if(b==2){
			offset_x = 0.5f;
			offset_y = 0f;
		}
		// offset to cell center
		else {
			offset_x = 0.5f;
			offset_y = 0.5f;
		}
		
		if (rkSteps==0){
		//**********************************************************************************
		// McCormack
			float phi_n_x, 		phi_n_y;
			float phi_n_hat_x, 	phi_n_hat_y;
			float phi_n1_hat_x,	phi_n1_hat_y;
			
			for (int i=1; i<=ssx; i++){
				for (int j=1; j<=ssy; j++){
					// generate particle at position [x|y]
					phi_n_x = (i+offset_x) ;
					phi_n_y = (j+offset_y) ;
					
					// send particle back to phi_n1_hat
					phi_n1_hat_x = phi_n_x -  dt * interpolate(phi_n_x		, phi_n_y-0.5f	,u) ;
					phi_n1_hat_y = phi_n_y -  dt * interpolate(phi_n_x-0.5f	, phi_n_y		,v) ;
					
					// send particle forward from phi_n1_hat to phi_n_hat
					phi_n_hat_x = phi_n1_hat_x +  dt * interpolate(phi_n1_hat_x		, phi_n_y-0.5f 	,u) ;
					phi_n_hat_y = phi_n1_hat_y +  dt * interpolate(phi_n1_hat_x-0.5f, phi_n_y		,v) ;
					
					// calculate phi_n1
					x =  phi_n1_hat_x + 0.5f * (phi_n_x - phi_n_hat_x);
					y =  phi_n1_hat_y + 0.5f * (phi_n_y - phi_n_hat_y);
					
					int lim_x = (int)(phi_n1_hat_x - offset_x);
					int lim_y = (int)(phi_n1_hat_y - offset_y);
					
					// Limiter - set result to range of first eulerstep.
					//float test = (float)interpolate(50,20, aOld);
					float l_1 = (float)interpolate(lim_x + 0,lim_y +0, aOld);
					float l_2 = (float)interpolate(lim_x + 1,lim_y +0, aOld);
					float l_3 = (float)interpolate(lim_x + 0,lim_y +1, aOld);
					float l_4 = (float)interpolate(lim_x + 1,lim_y +1, aOld);
					
					float max = Math.max(Math.max(Math.max(l_1, l_2),l_3),l_4);
					float min = Math.min(Math.min(Math.min(l_1, l_2),l_3),l_4);
					
					
					
					// interpolate the value of old field
					if(b==14) {
						a[i][j] = (float)interpolatePer(x-offset_x, y-offset_y, aOld);
						a[i][j] = Math.min( Math.max(a[i][j], min), max);
					}
					else{
						a[i][j] = (float)interpolate(x-offset_x, y-offset_y, aOld);
						a[i][j] = Math.min( Math.max(a[i][j], min), max);
					}
					
				}
			}
		}
		else{
		//**********************************************************************************
		// Runge Kutta
		
			for (int i=1; i<=ssx; i++){
				for (int j=1; j<=ssy; j++){
					// generate particle at position [x|y]
					x = (i+offset_x) ;
					y = (j+offset_y) ;
					
					// send particle backward at current timestep
					for(int k=0; k<rkSteps; k++){
						x = x -  (float)(1f/rkSteps) * dt * interpolate(x, y-0.5f, u);
						y = y -  (float)(1f/rkSteps) * dt * interpolate(x-0.5f, y, v);
					}
					
					// interpolate the value of old field
					if(b==4) 
						a[i][j] = (float)interpolatePer(x-offset_x, y-offset_y, aOld);
					else
						a[i][j] = (float)interpolate(x-offset_x, y-offset_y, aOld);
				}
			}
		}
		setBounds(b,a);
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
		
		setBounds(1, vF_u);
		setBounds(2, vF_v);
		
	}
	
	
	public void buoyancy(float[][] f, float k){
		float vpt,avpt;
		
		
		for (int i=0; i<=ssx; i++){
			for (int j=0; j<=ssy; j++){
				
				
				avpt = (float) (absT[j] * ( Math.pow( (p0/absP[j]) , 0.286)))*( 1 + 0.61f * qv[i][j]); 

				
				// B = (vpt-avpt)/ avpt - (g* qc)
				// vpt and avpt in Kelvin
				// g = 9.81 m/s²
				// qc in g/kg 
				vpt = pt[i][j] * ( 1 + 0.61f * qv[i][j]);
				f[i][j] = k*( ( (vpt-avpt) / avpt ) - 9.81f * qc[i][j] );
				d[i][j] =f[i][j];
				
			}
		}
		setBounds(2,f);
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
				
				exner = (float) ( Math.pow( (absP[j]/p0),0.286f )  );
				
				//p = (float) Math.pow(10*(1-(alt*tlr/t0)),(9.81/(tlr*rd)));
				//
				//compute 	T = pt[i,j]/( (^p/p)^k  )     
				// 			with  ^p = 100kPa   k = ~0.286
				//			T = pt[i,j]/( (100/p)^0.286)
				T =  pt[i][j]*exner;
				// conversion from Kelvin to Celsius
				T -= 273.15f;
				//
				//compute 	qs= (380.16/p)exp((17.67*T)/(T+243.5))
				// with T in °C and P in Pa
				//
				// qs nicht als feld
				
				
				qs = (float) ( (380.16f / (absP[j]*1000) ) * Math.exp( (17.67f * T) / (T + 243.5f) ) );
				
				d_qv  = Math.min(qs - qv[i][j],qc[i][j]);
				
				qv[i][j] = qv[i][j] + d_qv;
				qc[i][j] = qc[i][j] - d_qv;
				
				
				if(j>3 && qc[i][j]>qc_max){	
					qc_max=qc[i][j];
					System.out.println(1/qc_max);
				}
				
				//if(i==30 && j ==4)System.out.println("qs="+qs+" dqv="+d_qv+" qc="+qc[i][j]+" qv="+qv[i][j]);
				
				
				
				// Update the potential temperature according to condesation
				// Due to condensation latent energy is released in form of heat. -> change in pot temp
				//update potential Temperature ( Thermodynamics Equation )
				//delta p = L/(cp*PI)*(delta C)
				// L  = constant for latent heat released    	2501 J/kg
				// cp = specific heat capacity of dry air 		1005 J/(kg K)
				// C  = condensation rate = condensation per evaporation
				//	  = - min(qvs-qv, qc)  siehe waterCont = d_qv
				// PI = Exner Function = T/pt = 
				//_______________L_____/___cp___*________PI___________*____________C_________________________________
				T += 273.15f;
				
				
				pt[i][j] += (lh / ( cp * exner )) * (-d_qv);
				
				
				
			}
		}
		setBounds(4, qv);
		setBounds(5, qc);
		setBounds(3, pt);
		
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
	public void swapD(){
		temp = d;	
		d=dOld;		
		dOld=temp;
		}
	/**Swaps fields pt and ptOld */
	public void swapPt(){temp = pt;	pt=ptOld;	ptOld=temp;}
	/**Swaps fields qc and qcOld */
	public void swapQc(){temp = qc;	qc=qcOld;	qcOld=temp;}
	/**Swaps fields qv and qvOld */
	public void swapQv(){temp = qv;	qv=qvOld;	qvOld=temp;}
	

	
	/**
	 * 
	 * @param b 0=central Neumann || 1=uVel || 2=vVel || 3=potTemp || 4=waterVapor || 5=cloudWater
	 * @param f field to set boundaries
	 */
public void setBounds (int b, float[][] f){
	
	
		
		// b=0 central data neuman boundary
	if (b==0){
		for (int i=1 ; i<=ssx ; i++ ) { 
			f[i][0]     = f[i][1]; 
			f[i][ssy+1] = f[i][ssy]; 
		}
		for (int i=1 ; i<=ssy ; i++ ) { 
			f[0][i]     = f[1][i];   
			f[ssx+1][i] = f[ssx][i];
			
		}
	}

		
		// b=1 u velocity
		// bottom 	= noslip
		// top 		= free slip 
		// sides 	= user defined wind
		if (b==1){
			for (int i=1 ; i<=ssy ; i++ ) { 
				//f[1][i]    = -f[2][i];
				f[0][i]    = wind;   
				//f[1][i]    = wind; not sure wherer in or out
				f[ssx+1][i] = wind; 
			}
			for (int i=1 ; i<=ssx ; i++ ) { 
				f[i][0]   =  0; 
				f[i][ssy+1] =  f[i][ssy]; 
				
			}
		}
		
		// b=2 v
		// bottom 	= noslip
		// top 		= free slip 
		// sides 	= zero
		if (b==2){
			for (int i=1 ; i<=ssy ; i++ ) { 
				f[0][i]   = 0; 
				//f[1][i]   = 0; //not sure wherer in or out
			}
			for (int i=1 ; i<=ssx ; i++ ) { 
				f[i][0]     = 0; 
				//f[i][1]     = 0; 
				f[i][ssy+1] = 0; 
				 
			}
		}
		
		// b=3 potential temp
		// set to initial values
		// bottom noise
		if (b==3){
			float pt0=(float) (absT[0] * ( Math.pow( (p0/absP[0]) , 0.286)));  
			
			for(int i= 0; i<ssx+2; i++){
				f[i][0] = pt0;
				if(i>50 && i<150)
					f[i][2] = pt0+	5;//10-PerlinNoise.perlinNoise(i, time*0.8f+5000, 0.51f, 10f, 1f)*5f; 
				f[i][ssy+1] = 	(float) (absT[ssy+1] * Math.pow( (100/absP[ssy+1]) , 0.286 ) );  
			}
			for(int j= 0; j<ssy+2; j++){
				f[0][j]     = 	(float) (absT[j] * Math.pow( (100/absP[j] ) , 0.286));  
				f[ssx+1][j] = 	(float) (absT[j] * Math.pow( (100/absP[j] ) , 0.286));  
			}
			
		}
		
		// b=4 water vapor
		// periodic sides
		// top = 0
		// bottom = noise
		if (b==4){
			float qv0= (float) (hum	*(  (380/(absP[0]*1000)  ) * Math.exp( (17.67*(absT[0]-273.15)) / (absT[0]-273.15+243.5)))) ;
			float qv1= (float) (hum	*(  (380/(absP[ssy+1]*1000)  ) * Math.exp( (17.67*(absT[ssy+1]-273.15)) / (absT[ssy+1]-273.15+243.5)))) ;
			
			for(int i= 0; i<ssx+2; i++){
				f[i][ssy+1] = 	qv1;  
				f[i][0] 	= 	qv0;
				
				//PerlinNoise.perlinNoise(i, 50000000+time*0.8f, 0.51f, 10f, 1f)*0.9f; 
				if(i>50 && i<150)	f[i][0] =0.2f;//+= 	0.00041f;
			}
		}
		// b=5 cloud water
		// all to 0
		if (b==5){
			for (int i=0 ; i<ssy+2 ; i++ ) { 
				f[0][i]     = 0;
				f[ssx+1][i] = 0;
			}
			for (int i=0 ; i<ssx+2 ; i++ ) { 
				f[i][0]   	= 0; 
				f[i][ssy+1] = 0; 
			}
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
		
		//Boarder Conditions - not necessary!?
		if (xpos<0.5){xpos=0.5f;}
		if (xpos>ssx+0.5){xpos=ssx+0.5f;}
		if (ypos<0.5){ypos=0.5f;}
		if (ypos>ssy+0.5){ypos=ssy+0.5f;}
		//if (xpos<1){ypos=1f;}
		//if (xpos>=ssx+1){ypos=ssx+0.999999999999f;}
		//if (ypos<1){ypos=1f;}
		//if (ypos>=ssy+1){ypos=ssy+0.999999999999f;}
		
		// Sample positions
		x1 = (int) xpos;
		x2 = x1+1;
		y1 = (int) ypos;
		y2 = y1+1;
		
		// Distances
		dx = xpos%1;
		dy = ypos%1;
		//System.out.println(x1+" "+x2+" "+y1+" "+y2);
		// Interpolated Value
		return  	(1-dx) *( 	f[x1][y1]*(1-dy) 	+ 	f[x1][y2]*dy  	)
					+  dx  *(	f[x2][y1]*(1-dy)	+  	f[x2][y2]*dy  	);
		

		
	}
	
	
	public float interpolatePer(float xpos, float ypos, float[][] f){
		
		// Position Particle over boarder
		if(xpos<1){		xpos = (ssx+1) - (1-xpos);}
		else if(xpos>=((ssx+1))){	xpos = 1 + ( xpos - (ssx+1) );}
		
		// Sample positions
		x1 = (int) Math.floor(xpos);		
		x2 = x1+1;		
		
		//System.out.println("x "+x1);
		//System.out.println("x2 "+x2);
		
		//if particle is on boarder cell (between ssx and ssx+1)
		if (xpos>ssx){
			//System.out.println(xpos+" "+ypos);
		
			x1 = ssx;
			x2 = 1;
		}
		//System.out.println(xpos);
		if (ypos<0){ypos=0;} 
		if (ypos>=ssy+1){
			ypos=(ssy+0.99998f);}
		
		y1 = (int) Math.floor(ypos);
		y2 = y1+1;
		
		// Distances
		dx = (float) (xpos-Math.floor(xpos));//xpos%1;
		dy = (float) (ypos-Math.floor(ypos));//ypos%1;
		
		
		// Interpolated Value
		return  	(1-dx) *( 	f[x1][y1]*(1-dy) 	+ 	f[x1][y2]*dy  	)
					+  dx  *(	f[x2][y1]*(1-dy)	+  	f[x2][y2]*dy  	);
		

		
	}
	

	

	
}
