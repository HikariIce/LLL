import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.math.BigInteger;

public class GEN_KEY_PAIR {
	public void init() throws Exception{
		for(int i=0; i<1000; ++i){
			KeyPairGenrator kpg = new KeyPairGenrator();
			try{
				FileOutputStream out = new FileOutputStream("C:/Users/sun/Desktop/data/"+i+".data");
				FileOutputStream out1 = new FileOutputStream("C:/Users/sun/Desktop/data/"+i+".result");
				PrintStream p=new PrintStream(out);
				PrintStream p1=new PrintStream(out1);
	            p.println("[");
	            p.println("[1 0 "+kpg.e+" "+kpg.e+"]");
	            p.println("[0 1 "+kpg.n2+" "+kpg.n3+"]");
	            p.println("[0 0 0 0]");
	            p.println("[0 0 0 0]");
	            p.print("]");
	            p1.print(kpg.d);
			}catch(Exception e){
				e.printStackTrace();
			}
		}
//		saveToFile("public.key",kpg.n, kpg.e);
//		saveToFile("private.key", kpg.n, kpg.d);
	}
//	private void saveToFile(String fileName, BigInteger mod, BigInteger exp) throws IOException {
//		  ObjectOutputStream oout = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(fileName)));
//		  try {
//			  System.out.println("MOd  :"+mod);
//			  System.out.println("EXP  :"+exp);
//		    oout.writeObject(mod);
//		    oout.writeObject(exp);
//		  } catch (Exception e) {
//		    throw new IOException("Unexpected error", e);
//		  } finally {
//		    oout.close();
//		  }
//	}
	public static void main(String[] args) {
		try {
			new GEN_KEY_PAIR().init();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}