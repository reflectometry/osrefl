from pylab import cm, subplot, imshow, rc, arange, ones, figure, connect, subplots_adjust, axis, title, get_cmap, show, draw, close, ioff, ion
from numpy import outer


def show_colormap():
	ioff()
	rc('text', usetex=False)
	a=outer(arange(0,1,0.01),ones(10))
	fig = figure(figsize=(10,5))
	subplots_adjust(top=0.8,bottom=0.05,left=0.01,right=0.99)
	maps=[m for m in cm.datad.keys() if not m.endswith("_r")]
	maps.sort()
	l=len(maps)+1
	i=1

	for m in maps:
		subplot(1,l,i)
		axis("off")
		imshow(a,aspect='auto',cmap=get_cmap(m),origin="lower")
		title(m,rotation=90,fontsize=10)
		i=i+1
		
	show()
	draw()
	ion()
	print cname
	return cname

def change_colormap(image):
	''' takes matplotlib.image.AxesImage instance as argument, then
	displays a plot of all colormaps in cm module.  Click on one to 
	change the image colormap to the selection '''
	
	ioff()
	rc('text', usetex=False)
	a=outer(arange(0,1,0.01),ones(10))
	fig = figure(figsize=(10,5))
	subplots_adjust(top=0.8,bottom=0.05,left=0.01,right=0.99)
	maps=[m for m in cm.datad.keys() if not m.endswith("_r")]
	maps.sort()
	l=len(maps)+1
	i=1

	for m in maps:
		subplot(1,l,i)
		axis("off")
		imshow(a,aspect='auto',cmap=get_cmap(m),origin="lower")
		title(m,rotation=90,fontsize=10)
		i=i+1
		
	def clicker(event):
		if event.inaxes:
			cname = event.inaxes.title.get_text()
			print cname
			image.set_cmap(get_cmap(cname))
			fignum = image.figure.number
			fig = figure(fignum)
			show()
			draw()
		return
	
	connect('button_press_event',clicker)
	
	show()
	draw()
	ion()
	return

if __name__ == '__main__': show_colormap()