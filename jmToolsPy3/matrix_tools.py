"""
matrix_tools: Convenience functions for Programming the Matrix
==============================================================
	
 Ver	Date		Who	Comments
======	==========	===	================================================
0.0.90	2015-08-10	JRM	Initial prototype: complexPolarPlot
"""
# -*- coding: utf-8 -*-
def complexPolarPlot(a, color='#CC0000'):
	"""complexPolarPlot(a)
	Plot a list or vector of complex points as a polar plot.
	Defaults to red points.

	Input:
	a - a vector or list of points
	color - default red
	
	Returns:
	the fig"""
	import matplotlib.pyplot as plt
	import numpy as np
	for x in a:
		plt.polar(np.angle(x),np.abs(x), marker='.',color=color)
	# fig = plt.gcf()
	# return fig
