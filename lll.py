import os

if __name__ == '__main__':
	n = 0;
	dic = {}
	for i in xrange(1000):
		print str(i)+':'
		s = './lllv2 4 3 4 <./data/' + str(i) + '.data'
		val = os.popen(s).readline()[:-1]
		print 'LLL result : ',val
		if val[0] == '-':
			val = val[1:]

		with open('./data/' + str(i) + '.result','r') as f:
			res = f.readline()
		print 'RSA data d : ',res
		print 'Compare    : ',val==res
		if(val == res):
			n = n+1
		else:
			dic[i] = 'wrong --- ' + val + '\nright --- ' + res
		print '\n'
	
	print 'True : ',n
	print 'False: ',len(dic)
'''
	for (key,value) in dic.items():
		print key
		print value
		print '\n'	
'''	
	
	


	
