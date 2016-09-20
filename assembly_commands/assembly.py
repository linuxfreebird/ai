sequence = '3333333333322222222122222122222222222212222222222221222222222222222112222222222222222222222212222222222222221222222222222222222122222222222212222' #set of commands to execute in sequence to spell 'hello world'
listchr = 'abcdefghijklmnopqrstuvwxyz ' # allowed characters in a string to increment 
sinput = '' # starting sentence to alter
output = sinput + "" # altered sentence
pntmax = len(output)-1 # length of sentence
pnt = 0 # pointer to individual characters
for chr in sequence:
    if chr ==   '1':
        # increment pointer
        if pnt < pntmax:
            pnt += 1
        else:
            pnt = 0
    elif chr == '2':
        # increment char
        index = listchr.index(output[pnt])
        if index < len(listchr)-1:
            index += 1
        else:
            index = 0
        listoutput = list(output)
        listoutput[pnt] = listchr[index]
        output = "".join(listoutput)
    elif chr == '3':
        # add char
        output += " "
        pntmax += 1
    else:    # '4'
        # remove char
        if len(output) > 0:
            output = output[:-1]
            pntmax -= 1
print sinput
print output
