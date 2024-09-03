import argparse
import gzip

parser = argparse.ArgumentParser(description="hello")
parser.add_argument("r1", help="",type=str)
parser.add_argument("r2", help="", type=str)
parser.add_argument("res1", help="", type=str)
parser.add_argument("res2", help="", type=str)
parser.add_argument("r1Out", help="",type=str)
parser.add_argument("r2Out", help="", type=str)

args = parser.parse_args()




if args.r1.endswith(".gz"):
	r1Reader=gzip.open(args.r1)
else:
	r1Reader=open(args.r1)
if args.r2.endswith(".gz"):
	r2Reader=gzip.open(args.r2)
else:
	r2Reader=open(args.r2)


r1Writer=open(args.r1Out, "w")
r2Writer=open(args.r2Out, "w")

write=False

oldR1=""
oldR2=""

seqCount=0
lineCount=0
seqCountWritten=0

for line in r1Reader:
	line2 = r2Reader.readline()

	lineCount+=1

	if lineCount==1:
		oldR1=line
		oldR2=line2
	elif lineCount==2:
		if line.upper().startswith(args.res1) and line2.upper().startswith(args.res2):
			write=True
			seqCountWritten+=1
			r1Writer.write(oldR1)
			r2Writer.write(oldR2)
		#else:
		#	print(line[0:5])
		#	print(line2[0:5])

	if write==True:
		r1Writer.write(line)
		r2Writer.write(line2)


	if lineCount==4:
		seqCount+=1
		lineCount=0
		write=False
		oldR1=""
		oldR2=""


r1Reader.close()
r2Reader.close()
r1Writer.close()
r2Writer.close()


print("%s/%s sequences had correct restriction sites and were written to the new file!"%(seqCountWritten, seqCount))
