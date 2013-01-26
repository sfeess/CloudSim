import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;


public class Field {
	
	// linearize Array Position
	public static int flin(int i, int j, int ssx){
		return ((i)+(ssx+2)*(j));
	}
	
	// Circular field v
	public static float[] circFieldV(int x, int y, float t){
		float[] field = new float[(x+2)*(y+2)] ;
		for(int i=0; i<x+2; i++){
			for(int j=0; j<y+2; j++){ 
				//field[flin(i,j,x)]= (float) -Math.cos((Math.PI/(x))*i);
				field[flin(i,j,x)]=(float) (i-(x/2))/x;
			}
		}
		return field;
	}
	// Circular field u
	public static float[] circFieldU(int x, int y, float t){
		float[] field = new float[(x+2)*(y+2)] ;
		for(int i=0; i<x+2; i++){
			for(int j=0; j<y+2; j++){ 
				//field[flin(i,j,x)]= (float) (Math.cos((Math.PI/(y))*j));
				field[flin(i,j,x)]= (float) -(j-(y/2))/y;
			}
		}
		return field;
	}
	
	// Testfeld erstellen mit der größe x+2 y+2 - Verlauf quer
	public static float[] gradField(int x, int y){
		float[] field = new float[(x+2)*(y+2)] ;
		// Feld füllen
		for(int i=0; i<x+2; i++){
			for(int j=0; j<y+2; j++){
				//field[lin(i,j,x+2)] =(float)i/(float)(2*(x+2)) + (float)j/(float)(2*(y+2));
				field[flin(i,j,x)] = (float)i/(float)(x+2);
				//if(i>(x/2-x/6) && i<(x/2+x/6) && j>(y/2-y/6) && j<(y/2+y/6)){
				//field[flin(i,j,x)] = 0.6F;}//(float)(i%50)/(float)(50);}
			
			}
		}
		//for (int m=0; m<200; m+=10)
		//System.out.println("x:"+m+"y:"+m+" value:"+field[flin(m,m,x)]);

		return field;	
	}
	
	public static float[] boxField(int x, int y, float v){
		float[] field = new float[(x+2)*(y+2)] ;
		// Feld füllen
		for(int i=0; i<x+2; i++){
			for(int j=0; j<y+2; j++){

				if(i>(x/2-x/6) && i<(x/2+x/6) && j>(y/2-y/6) && j<(y/2+y/6)){
				field[flin(i,j,x)] = v;}
				else {field[flin(i,j,x)] =0F;}
			
			}
		}
		return field;	
	}
	
	public static float[] smlBoxField(int x, int y, float v){
		float[] field = new float[(x+2)*(y+2)] ;
		// Feld füllen
		for(int i=0; i<x+2; i++){
			for(int j=0; j<y+2; j++){

				if(i>(x/2-x/20) && i<(x/2+x/20) && j>(y*1/6-y/20) && j<(y*1/6+y/20)){
				field[flin(i,j,x)] = v;}
				else {field[flin(i,j,x)] =0F;}
			
			}
		}
		return field;	
	}
	
	public static float[] lineField(int x, int y){
		float[] field = new float[(x+2)*(y+2)] ;
		// Feld füllen
		
		for(int i=0; i<x+2; i++){
			for(int j=0; j<y+2; j++){

				if(i>(x/2-x/8) && i<(x/2+x/8)){
				field[flin(i,j,x)] = 1F;}
				else {field[flin(i,j,x)] =0F;
				
				if(i==0 || i==x+1 || j==0 || j==y+1)
				field[flin(i,j,x)] =0F;
				}
			
			}
		}
		return field;	
	}	
	
	public static float[] constField(int x, int y, float value){
		float[] field = new float[(x+2)*(y+2)] ;
		
		for(int i=0; i<x+2; i++){
			for(int j=0; j<y+2; j++){
				field[flin(i,j,x)]= 0;
			}
		}
		
		
		
		for(int i=1; i<x+1; i++){
			for(int j=1; j<y+1; j++){
				field[flin(i,j,x)]= value;
			}
		}
		return field;
		
	}
	
	public static float[] numbField(int x, int y){
	float[] field = new float[(x+2)*(y+2)] ;
		for(int i=0; i<x+2; i++){
			for(int j=0; j<y+2; j++){
				field[flin(i,j,x)]= (float)j/8f+(float)i/8f;
			}
		}
		return field;
	}

	public static float[] imgField(int x, int y){
		BufferedImage img = null;
		
		// Bild einlesen
		try {
			img = ImageIO.read(new File("smile.gif"));
		} catch (IOException e) {
			System.out.println("fehler");
			e.printStackTrace();
		}
		int h = img.getHeight(null);
		int w = img.getWidth(null);

		// Testfeld für bilddaten  anlegen
		float[] field = new float[(x+2)*(y+2)];
		WritableRaster raster= img.getRaster();
		
		// Bilddaten in testfeld schrieben
		for (int i=0; i<x+2; i++){
			for (int j=0; j<y+2; j++){
				int c = 0;
				if(i<w && j<h){
				c = raster.getSample(i,j, 0);}
				//c = raster.getSample(i,j, 0)-127;}
				//System.out.println(i+" "+j);
				field[flin(i,j,x)]=(float)c/(float)255;
			}
		}
		
		

		
		
		
		
		
		//float[] pixels = fs.evaluate(8,8,testfield);
		
		
		return field;
		}
	/**
	 * 
	 * @param x
	 * @param y
	 * @param p
	 * @param scale
	 * @param time
	 * @param amount
	 * @return
	 */
	public static float[] noiseEmitField(int x, int y, float p, float scale, float time, float amount){
		float[] field = new float[(x+2)*(y+2)] ;
		// Feld füllen
		
		for(int i=0; i<x+2; i++){
			for(int j=0; j<y+2; j++){

				
				if(j<2)field[flin(i,j,x)] = PerlinNoise.perlinNoise(i, j+time*0.2f, p, scale, 2f)*amount;
				
				if(i==0 || i==x+1 || j==0 || j==y+1)
				field[flin(i,j,x)] =0F;
				}
			
			}
		
		return field;	
	}	


}
