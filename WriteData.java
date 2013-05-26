import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import javax.imageio.ImageIO;


public class WriteData {
	
	public static void boxValuesOut(int x1, int x2, int y1, int y2){
		for (int i=x1; i<=x2; i++)
			for (int j=y1; j<=y2; j++)
				valuesOut(i, j);
		
	}

	
	public static void valuesOut(int x, int y){
		new File("values").mkdirs();
		
		File file = new File("values\\FluidValues_x"+x+"_y"+y+".txt");
		
	     try {
	       FileWriter writer = new FileWriter(file ,true);
	       
	       // Text wird in den Stream geschrieben
	       String t;
	       t = String.valueOf(FluidViewer.fs.step);
	       writer.write(t+",");
	       t = String.valueOf(FluidViewer.fs.u[x][y]);
	       writer.write(t+",");
	       t = String.valueOf(FluidViewer.fs.v[x][y]);
	       writer.write(t+",");
	       t = String.valueOf(FluidViewer.fs.pt[x][y]);
	       writer.write(t+",");
	       t = String.valueOf(FluidViewer.fs.qv[x][y]);
	       writer.write(t+",");
	       t = String.valueOf(FluidViewer.fs.qc[x][y]);
	       writer.write(t+",");
	       writer.write(System.getProperty("line.separator"));

	       writer.flush();
	       writer.close();
	    } 
	    catch (IOException e) {
	      e.printStackTrace();
	      System.out.println("Error");
	    }
	 }
	

	
	
	

	public static void imgOut(){
		
		new File("Cloud_Out").mkdirs();
		
		try {
			String x = String.valueOf(FluidViewer.fs.step);
			if (x.length()==1)x= "000"+x;
			if (x.length()==2)x= "00"+x;
			if (x.length()==3)x= "0"+x;
			ImageIO.write(FluidPanel.img, "bmp", new File("Cloud_Out\\CloudImg_"+x+".bmp"));
		} 
		catch (IOException ex) {
			System.out.println("Error");
        }
    }
	


}
