import javax.swing.JFrame;

public class ClientTest{
	public static void main(String[] args){
		Client you;
		you = new Client(args[0]);
		you.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		you.startRunning();
	}
}
