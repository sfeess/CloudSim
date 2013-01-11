import java.awt.BorderLayout;
import java.awt.Graphics;

import javax.swing.*;



public class Test2  extends JPanel {
	 
	   /**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public Test2() {              // Construktor der Klasse
	   }
	 
	   @Override
	   public void paintComponent(Graphics g) {
	     super.paintComponent(g);            
	     // meine Zeichnung
	     g.fillRect(100,100,150,250);
	   }
	 }