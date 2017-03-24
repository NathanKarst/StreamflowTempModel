import sys
import os
import numpy as np 


def main():
	for i in range(1,366):
		num = str(i)
		if len(num)==1:
			num = '00'+num
		elif len(num)==2:
			num = '0'+num
		else: 
			num = num

		a = 'diff_' + num + '.tif'
		b = 'beam_' + num + '.tif'
		command = "gdal_calc.py -A %s -B %s --outfile=rs_%s.tif --calc='A+B'"%(a,b,num)
		os.system(command)


if __name__ == '__main__':
	main()