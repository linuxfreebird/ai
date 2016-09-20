import sys
print "recursion limit = " + str(sys.getrecursionlimit())
from cmd import Cmd
remembered_line = ""
remembered_list = []
listchar ='0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ#!"$%&\'()*+,-./:;<=>?@[\\]^_`{|}~ \t\n\r\x0b\x0c'
totalchar = len(listchar)
network = [[for x in range(totalchar)] for y in range(totalchar)]
depth = 0
for x in network:
    for y in x:
        y = 0
def sorted_indices(s):
    return sorted(range(lens(s)),key=lambda k: s[k])
def recurse(indices,branchs,depth):


class MyPrompt(Cmd):
    def do_q(self, args):
        remembered_line = str(args)
        remembered_list = list(remembered_line)
        indices = []
        for s in remembered_list:
            indices.append(listchar.index(s))
        for i in range(len(indices)-1):
            network[indices[i]][indices[i+1]] += 1

    def do_quit(self, args):
        """Quits the program."""
        print "Quitting."
        raise SystemExit


if __name__ == '__main__':
    prompt = MyPrompt()
    prompt.prompt = '> '
    prompt.cmdloop('Starting prompt...')
