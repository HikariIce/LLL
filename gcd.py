import os

def	gcd(a,b):
	if b==0:
		return a
	return gcd(b,a%b)
	
def exgcd(a,b):
	if b==0:
		return (1,0)
	x,y = exgcd(b,a%b)
	x = x-a/b*y
	return (y,x)
	
def optimal(a,b,d,x,y):
	m = x*d%b
	m0 = m-b
	m2 = m+b
	dic1 = {}
	dic1[abs(m)] = m
	dic1[abs(m0)] = m0
	dic1[abs(m2)] = m2
	
	n = y*d%a
	n0 = n-a
	n2 = n+a
	dic2 = {}
	dic2[abs(n)] = n
	dic2[abs(n0)] = n0
	dic2[abs(n2)] = n2
	return dic1[min(dic1.keys())],dic2[min(dic2.keys())]

	
if __name__ == '__main__':
	n = 0;
	dic = {}
	
	for i in xrange(1000):
		s = './lllv2 4 3 4 <./data/' + str(i) + '.data'
		f = os.popen(s)
	
		a = f.readline()[:-1]
		if a[0] == '-':
			a = a[1:]
	
		b = f.readline()[:-1]
		if b[0] == '-':
			b = b[1:]
		
		a = int(a)
		b = int(b)
		
		with open('./data/' + str(i) + '.result','r') as f:
			d = f.readline()
		d = int(d)
		
		x,y = exgcd(a,b)
		x,y = optimal(a,b,d,x,y)
		print i,":",x,y
		dic[i] = (x,y)
	
	
'''
	for (key,value) in dic.items():
		print key
		print value
		print '\n'	
'''	
	
	
'''	
		print 'a : ',a
		print 'b : ',b
		if gcd(a,b)==1:
			n = n + 1
		print '\n'

	print n
'''
	
	

	
