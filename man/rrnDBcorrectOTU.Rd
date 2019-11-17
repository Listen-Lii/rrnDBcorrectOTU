\name{rco}
\alias{rco}
\title{use rrnDB database to coerrect 16S copy numbers for OTU table}
\usage{
rco(otu,classifer,rrnDB)
}
\arguments{
\item{otu}{OTU table}
\item{classifer}{RDP classifer results}
\item{rrnDB}{rrnDB database in RDP format}
}
\description{
input otu table,RDP classifer file, and rrnDB database (rdp format) in this order,we can get two outputs.
One is OTU table with corresponding 16S copy numbers information;
Another is corrected OTU tables.
}

