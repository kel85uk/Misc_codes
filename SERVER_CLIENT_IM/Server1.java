import java.io.*;
import java.net.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class Server1 extends JFrame{
  private JTextField userText;
  private JTextArea chatWindow;
  private ObjectOutputStream output;
  private ObjectInputStream input;
  private ServerSocket server;
  private Socket connection;

  // constructor
  public Server1(){
    super("Kelvin's Server");
    userText = new JTextField();
    userText.setEditable(false);
    userText.addActionListener(
      new ActionListener(){
          public void actionPerformed(ActionEvent event){
            sendMessage(event.getActionCommand());
            userText.setText("");
          }
      }
    );
    add(userText,BorderLayout.NORTH);
    chatWindow = new JTextArea();
    add(new JScrollPane(chatWindow));
    setSize(300,150);
    setVisible(true);
  }

  //set up and run server
  public void startRunning(){
    try{
      server = new ServerSocket(8080,100);
      while(true){
    	  try{
    		  //connect and have conversation
    		  waitForConnection();
    		  setupStreams();
    		  whileChatting();
    	  }catch(EOFException eofException){
    		  showMessage("\n Server ended the connection!");
    	  }finally{
    		  closeCrap();
    	  }
      }
    	      	  
    }catch(IOException ioException){
        ioException.printStackTrace();
    }
  }
  
  // wait for connection, then display connection information
  public void waitForConnection() throws IOException{
	  showMessage("Waiting for someone to connect... \n");
	  connection = server.accept();
	  showMessage("Now connected to " + connection.getInetAddress().getHostName() );
  }
  
  // set up streams
  private void setupStreams() throws IOException{
	  output = new ObjectOutputStream(connection.getOutputStream());
	  output.flush();
	  input = new ObjectInputStream(connection.getInputStream());
	  showMessage("\n Streams set up \n");
  }
  
  // Fun chatting
  private void whileChatting() throws IOException{
	  String message = " You are now connected! ";
	  sendMessage(message);
	  ableToType(true);
	  do{
		  //have conversation
		  try{
			  message = (String) input.readObject();
			  showMessage("\n" + message);
		  }catch(ClassNotFoundException classNotFoundException){
			  showMessage("\n idk wtf that user sent! \n");
		  }
	  }while(!message.equals("LEEPING: END"));
  }
  
  // Close crap after chatting
  private void closeCrap(){
	  showMessage("\n Closing connection... \n ");
	  ableToType(false);
	  try{
		  output.close();
		  input.close();
		  connection.close();
	  }catch(IOException ioException){
		  ioException.printStackTrace();
	  }
  }
  
  // Utility methods
  // Send message to client
  private void sendMessage(String message){
	  try{
		  output.writeObject("KELVIN: " + message);
		  output.flush();
		  showMessage("\nKELVIN: " + message);
	  }catch(IOException ioException){
		  chatWindow.append("\n ERROR: CANNOT SEND MESSAGE! \n");
	  }
  }
  
  // Update chat window
  private void showMessage(final String text){
	  SwingUtilities.invokeLater(
		  new Runnable(){
			  public void run(){
				  chatWindow.append(text);
			  }
		  }
	  );
  }
  
  // Able to type stuff
  private void ableToType(final boolean tof){
	  SwingUtilities.invokeLater(
		  new Runnable(){
			  public void run(){
				  userText.setEditable(tof);
			  }
		  }
	  );
  }
}
