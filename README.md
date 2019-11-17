# rrnDBcorrectOTU
cororect 16S rRNA copy numbers by rrnDB database for a OTU table

There are three input files:
OTU table
RDP classifer result file 
rrnDB database (RDP format) 
Input should be put in this order.
We can get two outputs. One is OTU table with corresponding 16S copy numbers information (whole.res); Another is corrected OTU tables(correct.table).

Usage
output = rco(otu,classifer,rrnDB)
output$cholw.res
output$correct.table

Attention!!!!
"-1" in the forth column (CopyNumber) in whole.res file means this taxa is not collected in rrnDB database.


