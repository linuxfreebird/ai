import subprocess
import string
import random
import os

listchar ='0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!"$%&\'()*+,-./:;<=>?@[\\]^_`{|}~ \t\n\r\x0b\x0c'

scount = 0
fcount = os.path.getsize('ai.py') - 1

while(True):
	try:
		subprocess.check_call(["python", "ai.py"])
		scount = fcount
	except subprocess.CalledProcessError:
		print 'fail'
	f = open('ai.py','r')
	txtf = f.read()
	listf=list(txtf)
	if (random.randint(0,4**(fcount-scount)) == 0 ):
		listf.append(random.choice(listchar))
		fcount = fcount + 1
	if(scount < fcount):
		listf[random.randint(scount,fcount)] = random.choice(listchar)
	f = open('ai.py','w')
	txtf = "".join(listf)
	f.write(txtf)
	f.close()
