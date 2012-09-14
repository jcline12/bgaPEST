
print 'I am the lizard king. I can do Anything!'

indat = open('test_data.dat','r').readlines()


ofp = open('test.out','w')
for line in reversed(indat):
	ofp.write(line)
ofp.close() 
