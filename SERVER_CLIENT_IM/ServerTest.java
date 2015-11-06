import javax.swing.JFrame;

public class ServerTest{
	public static void main(String[] args){
		Server1 me = new Server1();
		me.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		me.startRunning();
	}
}
