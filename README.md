The program is compiled by executing the make command.
The program reads a text file containing one IP direction per line (with no port). The file must be first converted to a binary format using the txtToBin program. The algorithm receives the binary file. 
The purpose of this is allowing the parallel reading of the file.

The algorithm presented here can be implemented in the p4 language in an architecture with pipelines. With 4 pipelines, the execution time is reduced by half.
