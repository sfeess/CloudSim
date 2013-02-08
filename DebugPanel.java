import java.awt.BorderLayout;
import java.awt.EventQueue;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import javax.swing.JLabel;
import java.awt.FlowLayout;
import javax.swing.BoxLayout;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import java.awt.GridLayout;
import javax.swing.JTextPane;
import javax.swing.JTable;
import java.awt.Color;
import java.awt.Window.Type;
import javax.swing.SwingConstants;


public class DebugPanel extends JPanel {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					DebugPanel frame = new DebugPanel();
					
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the frame.
	 */
	
	 JLabel lblDebugvalue1;
	 JLabel lblDebugvalue2;
	 JLabel lblDebugvalue3;
	 JLabel lblDebugvalue4;
	 JLabel lblDebugvalue5;
	 JLabel lblDebugvalue6;
	 JLabel lblDebugPos;
	
	public DebugPanel() {
		JPanel debugValues = new JPanel();
		debugValues.setBackground(Color.DARK_GRAY);
		debugValues.setSize(137,115);
		debugValues.setLayout(new GridLayout(0, 1, 0, 0));
		lblDebugPos = new JLabel("x:000 y:000");
		lblDebugPos.setForeground(Color.WHITE);
		debugValues.add(lblDebugPos);
		lblDebugvalue1 = new JLabel("v1");
		lblDebugvalue1.setHorizontalAlignment(SwingConstants.LEFT);
		lblDebugvalue1.setForeground(Color.WHITE);
		debugValues.add(lblDebugvalue1);
		lblDebugvalue2 = new JLabel("v2");
		lblDebugvalue2.setForeground(Color.WHITE);
		debugValues.add(lblDebugvalue2);
		lblDebugvalue3 = new JLabel("v3");
		lblDebugvalue3.setForeground(Color.WHITE);
		debugValues.add(lblDebugvalue3);
		lblDebugvalue4 = new JLabel("v4");
		lblDebugvalue4.setForeground(Color.WHITE);
		debugValues.add(lblDebugvalue4);
		lblDebugvalue5 = new JLabel("v5");
		lblDebugvalue5.setForeground(Color.WHITE);
		debugValues.add(lblDebugvalue5);
		lblDebugvalue6 = new JLabel("v6");
		lblDebugvalue6.setForeground(Color.WHITE);
		debugValues.add(lblDebugvalue6);
		
		
		
		
	}
	
	public void refresh(int x, int y){
		//x=40;		y=0;		
		lblDebugPos.setText("Position x: "+x+" y: "+y);
		lblDebugvalue1.setText(" u: "+FluidViewer.fs.u[x][y]);
		lblDebugvalue2.setText(" v: "+FluidViewer.fs.v[x][y]);
		lblDebugvalue3.setText(" pt: "+FluidViewer.fs.pt[x][y]);
		lblDebugvalue4.setText(" qc: "+FluidViewer.fs.qc[x][y]);
		lblDebugvalue5.setText(" qv: "+FluidViewer.fs.qv[x][y]);
		lblDebugvalue6.setText(" d: "+FluidViewer.fs.d[x][y]);
		
		
	}

}
