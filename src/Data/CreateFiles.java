package Data;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;
import java.net.*;
class CreateFiles{
	public static void main(String[] args){ 
		CreateFiles z = new CreateFiles();
		Create();
	}
	
public static void Create()
{
	BufferedWriter output = null;
	String filename=("Data");  
	File file = new File(filename);
	
    if(!file.exists()) 
	    file.mkdirs();
    
  //order search
    try {
    String order=("Data33_700000.txt");
    file = new File(order);
    BufferedWriter writer = new BufferedWriter(new FileWriter(file));
	FileWriter fw = new FileWriter(file,true); 
	
	for (int i=0;i<700000;i++)
	{
		//int i1=(int)(Math.random()*1000);       //  生成0-10000的随机数
		//int i2=(int)(Math.random()*1000);       //  生成0-10000的随机数
		String ss = i+" ";
		for ( int j = 0; j < 33; j ++)
			ss += (int)(Math.random()*1000) + " ";
		ss = ss.substring(0, ss.length()-1)+"\r\n";
		fw.write(ss);
		//int i2=Random.nextInt(100); 
		//fw.write(i +" "+i1+" " + i2 +"\r\n");
	}
	fw.close();
	writer.close();
	System.out.println("done");
	}
    catch ( IOException e ) {
        e.printStackTrace();
    } finally {
    }
	
	}
}
