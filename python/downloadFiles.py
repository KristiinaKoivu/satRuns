import urllib.request
import os
 
# files.txt contains the names of the files you want to download 
file1 = open('files.txt', 'r')
Lines = file1.readlines()

 
# here's the general path of the files to download 
url = 'https://a3s.fi/swift/v1/AUTH_70b34f161b3643938200c2ec96aa2ca0/PREBASruns/'
 
count = 0
# Strips the newline character
 
for line in Lines:

    urllib.request.urlretrieve((url+ urllib.parse.quote(line.strip())), os.path.basename(line).strip())
    #urllib.request.urlretrieve((url+ urllib.parse.quote(line.strip())), ("downloads/"+os.path.basename(line).strip()))
    print(os.path.basename(line).strip())
