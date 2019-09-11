

import matplotlib
from matplotlib.pylab import plt #load plot library
matplotlib.use('Agg')
from sys import stdout
import sys


filepath = str(sys.argv[1])
print('Reading file ', filepath)

# Frag,562480,166774007,562667,166774194,f,0,187,66,66,0.352941,0.352941,0,0
#SeqX length : 155270560
#SeqY length : 171031299


plt.figure(1, figsize = (8.5,11))
axes = plt.gca()
plt.gca().set_aspect('equal')

x_size = 1000
y_size = 1000

axes.set_xlim(0, x_size)
axes.set_ylim(0, y_size)

#ratio_x = x_size / 155270560
#ratio_y = y_size / 171031299


count = 0
with open(filepath) as openfileobject:
	for line in openfileobject:

		if(count == 6):
			split = [x.strip() for x in line.split(':')]
			seq_x_len = float(split[1])
			print("X len:", seq_x_len)

		if(count == 7):
			split = [x.strip() for x in line.split(':')]
			seq_y_len = float(split[1])
			print("Y len:", seq_y_len)

			ratio_x = x_size / seq_x_len
			ratio_y = y_size / seq_y_len
			print('Ratios: ', ratio_x, ratio_y)
		
		if((count+1) % 1000 == 0):
			stdout.write("\r%d" % count)
		if(count > 17):

			split = [x.strip() for x in line.split(',')]
			xdata = [float(split[1])*ratio_x, float(split[3])*ratio_x]
			ydata = [float(split[2])*ratio_y, float(split[4])*ratio_y]
			plt.plot(xdata, ydata, color="black")
			

		count = count + 1
		
plt.savefig('comparison.png')

print('\nSaved figure')



#plt.plot(x,y)

