import string

# Will expand this to list containing all lisp functions
lispfuncs = ['blah']

# Represents list of characters in recently produced code
uc = ['a', 'b', 'l', 'a', 'h', 'c', 'z']

# Eventually move this to loop and cycle though list of lisp functions
ws = list(lispfuncs[0])

fs = len(ws)

# Used as a boolean condition check, length is size of current lisp function we're searching for
bc = []
ebc = 1
for n in range(0, fs):
	bc.append(0)
	
for i in range(0, len(uc)):

	for k in range(0, fs):
		bc[k] = 0

	ebc = 1
	
	if i + fs < len(uc):
# Checking for identical match, mark each identical spot with 1		
		for j in range(i, i + fs):
			if uc[j] == ws[j - i]:
				print "%s    %s    %d" % (uc[j], ws[j-i], j)
				bc[j - i] = 1
		for m in range(0, fs):
			ebc = ebc*bc[m]
		
		if ebc == 1	:
			print "We have a match between element %d and %d" % (i, i+fs-1)
		
		print "Value of i is %d" % i