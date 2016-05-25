import java.math.BigInteger;


public class FindXY {
	public static void main(String[] args) {
		BigInteger res = new BigInteger("874192903529632456430885449811744750386939558049458917065474482895882917221941");
		BigInteger a   = new BigInteger("97082379893510002585565329776545598590802263690670536704258593656002576039284");
		BigInteger b   = new BigInteger("194616244275062438331928141375925560251323712214765160135664327303864884946953");

		BigInteger x = new BigInteger("7");
		BigInteger y = new BigInteger("1");
		System.out.println(x.multiply(a).add(y.multiply(b)).compareTo(res));
		
		System.out.println(gcd(a,b));
		System.out.println("x: "+exgcd(a,b)[0]+"  y: "+exgcd(a,b)[1]);
		System.out.println();
		
		
		BigInteger aa = exgcd(a,b)[0];
		BigInteger bb = exgcd(a,b)[1];
//		System.out.println(aa.multiply(a).add(bb.multiply(b)));
//		System.out.println();
		
		comp(b,res,aa);
		comp(a,res,bb);
	}
	
	public static BigInteger gcd(BigInteger a, BigInteger b){
		while(b.compareTo(BigInteger.ZERO)!=0){
			BigInteger r = b;
			b = a.divideAndRemainder(b)[1];
			a = r;
		}
		return a;
	}
	

	public static BigInteger[] exgcd(BigInteger a, BigInteger b){
		if(b.compareTo(BigInteger.ZERO)==0){
			BigInteger[] ret = new BigInteger[2];  //x : ret[0],   y : ret[1]
			ret[0] = BigInteger.ONE;
			ret[1] = BigInteger.ZERO;
			return ret;
		}
		BigInteger[] d = exgcd(b,a.divideAndRemainder(b)[1]);
		
		d[0] = d[0].subtract(a.divide(b).multiply(d[1]));
		
		BigInteger temp = d[0];
		d[0] = d[1];
		d[1] = temp;

		return d;
	}
	
	public static void comp(BigInteger b,BigInteger res,BigInteger x){
		BigInteger m = x.multiply(res).divideAndRemainder(b)[1];
		
		if(m.abs().compareTo(m.add(b).abs())<0)
			if(m.compareTo(m.subtract(b).abs())<0)
				System.out.println(m);
			else
				System.out.println(m.subtract(b));
		else
			if(m.add(b).compareTo(m.subtract(b).abs())<0)
				System.out.println(m.add(b));
			else
				System.out.println(m.subtract(b));
	}
	
	
}
