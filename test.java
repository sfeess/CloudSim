import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;


public class test {
	static FluidSolver fs = new FluidSolver();
	public static void main(String[] args) {
		float[] f = Field.gradField(8,8);
	
		FluidSolver fs= new FluidSolver();
		fs.interpolate(3.1F, 3.1F, f);
	
	
	}
	
	
}
