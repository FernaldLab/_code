makeWaves = function( increment=0.01,
					  amplitude=1,
					  period=1,
					  shift=0,
					  series=c(1,6/5,4/3,3/2,9/5,2),
					  func='sin'
					  ) {
	x = seq(0, period*2*pi, by=increment);
	#x = x-shift;
	y = amplitude*eval(call(func,x));
	plot(x, max(series)*y ,type='n');
	for (r in series) {
		lines(x, r*y, type='l');
		x=x+shift;
	}
	
}