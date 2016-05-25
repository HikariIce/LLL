import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.Random;


public class KeyPairGenrator {
	public BigInteger p,q,n,e,d,v;
	public BigDecimal nn,n2,n3;
	Random rnd,rnd1;
	int prefferd_block_size;
	public KeyPairGenrator() {
		rnd=new Random();
		rnd1=new Random();
		p=new BigInteger(512, 100,rnd );
		q=p.nextProbablePrime();
		n=q.multiply(p);
		v=(p.subtract(BigInteger.ONE).multiply(q.subtract(BigInteger.ONE)));
		d=getPublic(v);
		e=getPrivate(d,v);
		prefferd_block_size=(int)Math.min(Math.floor((1024-1)/8),256);
		
		System.out.println("n: " + n);
		System.out.println("v: " + v);
		System.out.println("p: " + p);
		System.out.println("q: " + q);
		System.out.println("e: " + e);
		System.out.println("d: " + d);
		
		System.out.println();
		nn = new BigDecimal(n);
		BigDecimal sqrtn = new BigDecimal(Math.sqrt(n.doubleValue())).setScale(0, BigDecimal.ROUND_HALF_UP);
		n2 = nn.subtract(sqrtn.multiply(new BigDecimal("1.5"))).setScale(0, BigDecimal.ROUND_HALF_UP);
		n3 = nn.subtract(sqrtn.multiply(new BigDecimal("2.5"))).setScale(0, BigDecimal.ROUND_HALF_UP);
		System.out.println(e.multiply(d).divideAndRemainder(v)[1]);
	}
	
	private BigInteger getPublic(BigInteger v2) {
		BigInteger ret =BigInteger.probablePrime((int)(v2.bitLength()/4)+4,rnd); //control the size of d.
		while(!(ret.gcd(v2).equals(BigInteger.ONE))){
		//System.out.println(ret);
		ret =ret.nextProbablePrime();}

		return ret;
	}
	
	
	
	
	public KeyPairGenrator(int keyLen) {
		rnd=new Random();
		rnd1=new Random();
		p=new BigInteger(keyLen/2, 100,rnd );
		q=new BigInteger(keyLen/2, 100,rnd1 );
		n=q.multiply(p);
		v=(p.subtract(BigInteger.ONE).multiply(q.subtract(BigInteger.ONE)));
		d=getPublic(v);
		e=getPrivate(d,v);
		prefferd_block_size=(int)Math.min(Math.floor((keyLen-1)/8),256);
	}
	private BigInteger getPrivate(BigInteger e2, BigInteger v2) {
		BigInteger i=BigInteger.valueOf(1);
		i=e2.modInverse(v2);
	
		return i;	
	}

	private BigInteger getPrivate1(BigInteger e2, BigInteger v2) {
		BigInteger i=BigInteger.valueOf(1);
		BigInteger toB1=(i.multiply(e2).mod(v2));
		while(!(toB1.equals(BigInteger.ONE))){
			i=i.add(BigInteger.ONE);
			toB1=i.multiply(e2);
			toB1=toB1.mod(v2);
			//System.out.println(i);
		}
		return i;
	}
	
	public BigInteger getPrivateKey(){return d;}
	public BigInteger getPublicKey(){return e;}
	
	public static void main(String[] args) {
		new KeyPairGenrator();
	}
}
