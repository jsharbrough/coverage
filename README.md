# coverage

This repository is designed to quickly obtain coverage levels of reads mapped to a reference 
sequence from a supplied vcf (Variant Call Format) file. If a region file is also provided, then the 
program will calculate mean and median coverage levels within the region(s) specified. If a region 
file is not specified, then the program will output mean and median coverage levels genome wide. 
Future iterations will include the ability to filter out high coverage regions, which should be 
especially helpful when determining average genome-wide coverage levels. The region file should be 
tab-delimited and formatted as follows:

Gene	Scaffold	Start	Stop

This region file is automatically generated and properly formated by the samPositions.py module 
from a sam file, but can also be generated manually. There should be no header in the regions file. 
Gene names are not required and can be identical; however, a tab must remain in order for the file 
to be parsed correctly. 

coverage.py will first attempt to open an index of the vcf file. If one does not exist in the 
current working directory, then the program will create one. Otherwise, the program will be able 
to load the stored index and quickly access coverage information stored in a vcf file. Output is to 
the standard out unless otherwise specified and will be tab-delimited.

USAGE

python coverage.py foo.vcf regions.txt > coverage.txt

Questions and comments should be addressed to Joel Sharbrough at jsharbro[at]gmail.com.

Copyright 2017 Joel Sharbrough under the MIT license. This program is provided as is and may be shared, edited, and used freely for non-commercial purposes ONLY.

All rights reserved.
