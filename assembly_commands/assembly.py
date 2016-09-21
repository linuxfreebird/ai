def assembly(sinput,sequence,listchr):
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
            output += output[0]
            pntmax += 1
        elif chr == '4':
            # remove char
            if len(output) > 0:
                output = output[:-1]
                pntmax -= 1
        else:
            0 + 0 # do nothing
    return output
print assembly('432331','123443243243241224312432','1234')
memory = [["" for y in range(8)] for x in range(8)]
print memory
