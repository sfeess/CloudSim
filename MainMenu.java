import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;


public class MainMenu{


}

class Exit implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			System.exit(0);
		}
	}

class Restart implements ActionListener {
	public void actionPerformed(ActionEvent e) {
		
		// FluidViewer.reset();
	}
}