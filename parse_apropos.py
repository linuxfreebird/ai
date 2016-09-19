def letters(input):
    return ''.join(filter(str.isalpha, input))
file = open("apropos2.text",'r')
text = file.read()
list1 = text.split("\n")
list2 = []

for x in list1:
	list2.append(x.split())
types = []
for x in list2:
	if(len(x) == 0):
		types.append('blank')
        elif (len(x)==1):
		types.append('null')
	else:
                temp = ""
                for y in x:
                        if letters(y).islower() and ("|" not in y) :
                                temp = temp + y + " "
                types.append(temp) 
unique_types = list(set(types))
unique_types.sort()
sorted_list =[[] for i in range(len(unique_types))]
for i,y in enumerate(unique_types):
        for j,x in enumerate(types):
                if x == y :
                        temp = ""
                        for z in list2[j]:
                                if not (letters(z).islower() and ("|" not in z)):
                                        temp = temp + z + " "
                        sorted_list[i].append(temp)
for i,title in enumerate(unique_types):
        with open('./apropos/' + title + '.txt','w') as f:
                f.write("\n".join(sorted_list[i]))
