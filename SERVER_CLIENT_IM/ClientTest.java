import javax.swing.JFrame;

public class ClientTest{
	public static void main(String[] args){
		Client you;
		you = new Client("127.0.0.1");
		you.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		you.startRunning();
	}
}