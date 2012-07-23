In order to compile, just type:

make

from terminal.

To execute the program, type:

./lfm

--------------------------------
The program will ask for the name of a file where it can read the list of edges and their weights (in fact the program works on weighted networks).
This file is expected like this:

source_node target_node weight


Actually, since the program works only on undirected graphs, it does not care about who is the source and who is the target. Repetitions and self-loops are ignored.

The node is supposed to be an integer number, and the weight a real one.


If you do not want to type from terminal, you could edit the program here:

./Library/lfm.cpp, line 55

Just below this line, there are some parameters that you might like to change as well...
(If you do not know c++, email me please)
-------------------------------



The program will write the results in a file called "output.dat". The several collections of overlapping nodes are ordered according to their relative frequency.





