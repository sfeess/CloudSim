
public class FS {
	
	//Gridsize (without +2 Boarder)
	int ssx, ssy;
	
	//density
	float[] d, dOld;
	
	//output
	float[] out;
	
	//velocity
	float[] u, uOld;
	float[] v, vOld;
	
	float time;
	int size;
	float dt;
	
	// evaluate liest pixelwert aus feld aus
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
				
				//transform Pixle positions in Cellspace 
				float xpos= xscale*i;
				float ypos= yscale*j;
				
				//System.out.println(xpos+" "+ypos);
				
				//Samplepixel x1=pos matching pixel  x2= pos neighbor left or right
				int x1 = (int) (xpos); 
				int x2 = (int) Math.abs(Math.floor(xpos)+(int)Math.signum(xpos%1-0.5));
				int y1 = (int) ypos;
				int y2 = (int) Math.abs(Math.floor(ypos)+(int)Math.signum(ypos%1-0.5));
				
				float dx = (float) Math.abs((xpos%1)-0.5);
				float dy = (float) (1-Math.abs((ypos%1)-0.5));
				
				//11= x1 y2
				//12= x2 y2
				//21= x1 y1
				//22= x2 y1

				//System.out.println(lin(i,j));
				//System.out.print(" ");
				//System.out.println(j);
				mapped[lin(i,j,(int)x)]=	
									f[lin(x1,y2)]*(1-dx)*(1-dy)+
									f[lin(x2,y2)]*dx*(1-dy)+
									f[lin(x1,y1)]*dy*(1-dx)+
									f[lin(x2,y1)]*dx*dy;
				
				
			}
			
		}
		return mapped;
	}
	
	
	// step führt solver zeitschrit aus
	public void step(){
		time+=dt;
		
		// Testfeldupdate (wieder löschen)
		u=Field.circFieldU(ssx,ssy,time);
		v=Field.circFieldV(ssx,ssy,time);
		out = Field.imgField(ssx,ssy);
		
		solveVel();
		solveDens();
	}
	
	//wandle 2D array pos in 1D array pos
	public int lin(int i, int j){
		return ((i)+(this.ssx+2)*(j));
	}
	public static int lin(int i, int j, int width){
		return ((i)+(width)*(j));
	}

	
	
	// solver setup
	public void setup(int sx, int sy, float dt){
		this.ssx = sx;
		this.ssy = sy;
		this.size = (sx+2)*(sy+2);
		this.dt = dt;
		
		
		
		reset();
		
	}
		
	//create zero-fields
	public void reset(){
		u = uOld = v = vOld = d = dOld = out = new float[size];
		for(int x= 0; x<size; x++){
			u[x] = uOld[x] = v[x] = vOld[x] = d[x] = dOld[x] = out[x] = 0.0f;
		}
	}
	
	//swap Fields
	public void swapU(){
		float[] temp = u;
		u=uOld;
		uOld=temp;
	}
	public void swapV(){
		float[] temp = v;
		v=vOld;
		vOld=temp;
	}
	public void swapD(){
		float[] temp = d;
		d=dOld;
		dOld=temp;
	}
	
	
	
	
	public void addForce(){}
	
	
	public void advect(){
		// advect()für jede dimension seperat?
		// advect (u an vel[u,v])
		// advect (v an vel[u,v])
		// advect (d an vel[u,v])
		
		// dOld[] & d[] = field to advect - Density or Velocity
		// uOld[] u[] vOld[] v[] = Velocity
		//
		
		// trace linear einen Timestep im raum
		// gesuchter Punkt =  aktueller Punkt - dt * velocity 
		// gesuchte xKoordinate x = aktuell PunktX - dt * velocityU
		// gesuchte yKoordinate y = aktuell PunktY - dt * velocityV
		// Suche Value im Feld [(int)x,(int)y]
		
		// interpoliere zu umliegenden Feldern
		// suche umliegenden Felder (wie evaluate) von pos[0,0] aus
		// betrachte nachkommastellen des neuen punktes und 
		// gewichte/interpoliere demnach den wert der umliegenden felder
		
	}
	
	
	public void diffuse(){}
	public void project(){}
	
	public void addSource(){}
	
	
	public void solveVel(){
		addForce();
		
		advect();
		diffuse();
		project();
		}
	
	public void solveDens(){
		
		addSource();
		diffuse();
		advect();
	}
	
	
	
	
	
	//
	
	

	
	
	 
	
}

























public class Test2{
	
	public static int lin(int i, int j, int width){
		return ((i)+(width)*(j));
	}
	static int x=2;
	static int y=5;
	
	public static void advect(int N, int b, float[] d, float[] d0, float[] u, float[] v, float dt){
		int i, j, i0, j0, i1, j1;  
		float x, y, s0, t0, s1, t1, dt0;  
		 
		dt0 = dt*N;  
		for ( i=1 ; i<=N ; i++ ) {   
			for ( j=1 ; j<=N ; j++ ) { 
				
				//vorherigen Koordinaten berechnen
				x = i-dt0*u[IX(i,j)]; 
				y = j-dt0*v[IX(i,j)];
				
				//bei boardern x und y vor grenzen setzen
				if (x<0.5) x=0.5F; 
				if (x>N+0.5) x=N+0.5F; 
				
				//Nachbarzellenpos in x Richtung
				i0=(int)x; 	i1=i0+1;    
				
				if (y<0.5) y=0.5F; 
				if (y>N+0.5) y=N+0.5F; 
				
				//Nachbarzellenpos in y Richtung
				j0=(int)y; j1=j0+1;    
				
				
				//s1 & t1 = nachkommastellen von x und y jeweils  
				s1 = x-i0; s0 = 1-s1; 
				t1 = y-j0; t0 = 1-t1;    
				
				d[IX(i,j)] = 	s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+         
								s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);   
				}  
			}  
		
		set_bnd ( N, b, d ); 
	
	} 
	}
	
	
	public static void main(String[] args){
		
		double x = Math.PI/2;
		System.out.println(Math.cos(x));
		
		
		
		
	}


	
}
