


public class PerlinNoise {
	
	
	


	// generate random numbers for x and y
	static public double noise1(int x, int y){
		int n = x + y * 57;
		n= (n<<13) ^ n;
		return ( 1.0 - ( (n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0);    
	}
	static public double noise2(int x, int y) {
		int n = x + y * 57;
		n= (n<<13) ^ n;
		return ( 1.0 - ( (n * (n * n * 15733 + 789227) + 1376312627) & 0x7fffffff) / 1073741827.0);    
	}
	static public double noise3(int x, int y) {
		int n = x + y * 57;
		n= (n<<13) ^ n;
		return ( 1.0 - ( (n * (n * n * 15737 + 789251) + 1376312629) & 0x7fffffff) / 1073741831.0);    
	}
	static public double noise4(int x, int y) {
		int n = x + y * 57;
		n= (n<<13) ^ n;
		return ( 1.0 - ( (n * (n * n * 15739 + 789311) + 1376312657) & 0x7fffffff) / 1073741833.0);    
	}

	// Smooth by interpolation with neighbor values
	public static float smoothN1(int x, int y)
	{
		float corners = (float) (( noise1(x-1, y-1)+noise1(x+1, y-1)+noise1(x-1, y+1)+noise1(x+1, y+1) ) / 16);
		float sides   = (float) (( noise1(x-1, y)  +noise1(x+1, y)  +noise1(x, y-1)  +noise1(x, y+1) ) /  8);
		float center  =  (float) (noise1(x, y) / 4);
		//return (corners + sides + center);
		return (float) noise1(x,y);
	}
	public static float smoothN2(int x, int y)
	{
		float corners = (float) (( noise2(x-1, y-1)+noise2(x+1, y-1)+noise2(x-1, y+1)+noise2(x+1, y+1) ) / 16);
		float sides   = (float) (( noise2(x-1, y)  +noise2(x+1, y)  +noise2(x, y-1)  +noise2(x, y+1) ) /  8);
		float center  =  (float) (noise2(x, y) / 4);
		return (corners + sides + center);
	}
	public static float smoothN3(int x, int y)
	{
		float corners = (float) (( noise3(x-1, y-1)+noise3(x+1, y-1)+noise3(x-1, y+1)+noise3(x+1, y+1) ) / 16);
		float sides   = (float) (( noise3(x-1, y)  +noise3(x+1, y)  +noise3(x, y-1)  +noise3(x, y+1) ) /  8);
		float center  =  (float) (noise3(x, y) / 4);
		return (corners + sides + center);
	}
	public static float smoothN4(int x, int y)
	{
		float corners = (float) (( noise4(x-1, y-1)+noise4(x+1, y-1)+noise4(x-1, y+1)+noise4(x+1, y+1) ) / 16);
		float sides   = (float) (( noise4(x-1, y)  +noise4(x+1, y)  +noise4(x, y-1)  +noise4(x, y+1) ) /  8);
		float center  =  (float) (noise4(x, y) / 4);
		return (corners + sides + center);
	}
	
	// Cosinus interpolate between a and b at point d     d=1 =>100% b   d=0 => 100% a
	static float interpolCos(float a, float b, float d) {
		float ft = d * 3.1415927f;
		float f = (float) ((1 - Math.cos(ft)) * .5);
		return  a*(1-f) + b*f;
	}
	
	// Interpolates between the smoothed Noisevalues
	public static float interpNoise1(float x, float y)
	{
		float frac_x = x - (int)x;
		float frac_y = y - (int)y;
		float v1 = smoothN1((int)x,     (int)y);
		float v2 = smoothN1((int)x + 1, (int)y);
		float v3 = smoothN1((int)x,     (int)y + 1);
		float v4 = smoothN1((int)x + 1, (int)y + 1);

		float i1 = interpolCos(v1 , v2 , frac_x);
		float i2 = interpolCos(v3 , v4 , frac_x);

		return interpolCos(i1 , i2 , frac_y);		 
	}
	
	
	
	public static float perlinNoise(float x, float y, float per,float scale, float contrast){
		// resulting value and persitence
		float res=0;
		// Number of octaves
		int oct = 5;
		
		res=0;
		
		for(int i=0; i<oct; i++){
			float freq = (float) Math.pow(2,i);
			float amp =  (float) Math.pow(per,i);
			res= res+interpNoise1(x*freq*1/scale, y*freq*1/scale)*amp;
		}
		
		res = res*contrast;
		
		res = Math.max(res,-1);
		res= Math.min(1, res);
		res = res+1;
		res = res*0.5f;
		res=Math.abs(res);
		
		return res;
	}
	
	
	public static void main(String args[]){
		
		for(int i=0; i<100;i++)	for(int j=0; j<100;j++) System.out.println(perlinNoise(i,j,0.7f,1,1));
	}

}
