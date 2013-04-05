'''
Created on Apr 4, 2013

@author: Ankit Sharma
'''
import re;

#str = "133149    165231,165260    2000,2000,2001,2001,2001,2001,2002,2002";
#words = str.split('\t',1);
#print(words[0]=='133149');
#print(words[0]==133149);

output_file = 'DBLP_Hadoop_Output.txt';
f = open(output_file,'w');

filename = 'DBLP_Hadoop_Input.txt';
lines = tuple(open(filename, 'r'));
#print(lines);
for line in lines:
    words = line.split('\t');
    years = words[2].split(',');
    authors = words[1].split(',');
    edge = words[0];
    
    if '-1' not in years[0]:       
        for author in authors:
            for year in years:
                f.write(edge.replace('\n','')+','+author.replace('\n','')+','+year.replace('\n','')+'\n');
                #print(edge.replace('\n','')+','+author.replace('\n','')+','+year.replace('\n','')+'');


f.close();