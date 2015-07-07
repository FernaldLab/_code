par(mfrow=c(2, 3));
labels = c('Day 0', 'Day 1', 'Day 2');
ylab = 'Aggressive beh / min';
col = 'grey';
pch = 21;
ylim = c(0, 8);
cex = 1.5; cex.axis = 1.5; cex.lab = 1.5;

###
toPlot = datNP[,11:13];
main = 'NP winners';

boxplot(toPlot, names = labels, ylab = ylab, main = main, ylim = ylim, cex.axis = cex.axis, cex.lab = cex.lab);
stripchart(toPlot, add = T, vertical = T, pch = pch, bg = col, cex = cex);
for (row in 1:nrow(toPlot))
{
	segments(1, toPlot[row, 1], 2, toPlot[row, 2], col = col);
	segments(2, toPlot[row, 2], 3, toPlot[row, 3], col = 'red');
}
rm(row, toPlot);

###
toPlot = datP[,11:13];
main = 'P winners';

boxplot(toPlot, names = labels, ylab = ylab, main = main, ylim = ylim, cex.axis = cex.axis, cex.lab = cex.lab);
stripchart(toPlot, add = T, vertical = T, pch = pch, bg = col, cex = cex);
for (row in 1:nrow(toPlot))
{
	segments(1, toPlot[row, 1], 2, toPlot[row, 2], col = col);
	segments(2, toPlot[row, 2], 3, toPlot[row, 3], col = 'red');
}
rm(row, toPlot);

###
toPlot = datMatched[,10:12];
main = 'All age matched';

boxplot(toPlot, names = labels, ylab = ylab, main = main, ylim = ylim, cex.axis = cex.axis, cex.lab = cex.lab);
stripchart(toPlot, add = T, vertical = T, pch = pch, bg = col, cex = cex);
for (row in 1:nrow(toPlot))
{
	segments(1, toPlot[row, 1], 2, toPlot[row, 2], col = col);
	segments(2, toPlot[row, 2], 3, toPlot[row, 3], col = 'red');
}
rm(row, toPlot);

###
toPlot = datNPnomatch[,11:13];
main = 'NP older';

boxplot(toPlot, names = labels, ylab = ylab, main = main, ylim = ylim, cex.axis = cex.axis, cex.lab = cex.lab);
stripchart(toPlot, add = T, vertical = T, pch = pch, bg = col, cex = cex);
for (row in 1:nrow(toPlot))
{
	segments(1, toPlot[row, 1], 2, toPlot[row, 2], col = col);
	segments(2, toPlot[row, 2], 3, toPlot[row, 3], col = 'red');
}
rm(row, toPlot);

###
toPlot = datNPmatch[,11:13];
main = 'NP age matched';

boxplot(toPlot, names = labels, ylab = ylab, main = main, ylim = ylim, cex.axis = cex.axis, cex.lab = cex.lab);
stripchart(toPlot, add = T, vertical = T, pch = pch, bg = col, cex = cex);
for (row in 1:nrow(toPlot))
{
	segments(1, toPlot[row, 1], 2, toPlot[row, 2], col = col);
	segments(2, toPlot[row, 2], 3, toPlot[row, 3], col = 'red');
}
rm(row, toPlot);

########################
###################
#############

plot(c(0, 1, 2, 
	   0, 1, 2,
	   0, 1, 2,
	   0, 1, 2,
	   0, 1, 2
	   ), 
	 c(apply(datNP[,11:13], 2, mean), 
	   apply(datP[,11:13], 2, mean),
	   apply(datMatched[,10:12], 2, mean),
	   apply(datNPnomatch[,11:13], 2, mean),
	   apply(datNPmatch[,11:13], 2, mean)
	   ),
	  ylab = 'Aggressive beh / min',
	  xlab = '',
	  axes = F,
	  type = 'n'
	 );

axis(at = c(0, 1, 2), side = 1, labels = c('Day 0', 'Day 1', 'Day 2'));
axis(side = 2);

lwd = 3;
lines(c(0, 1, 2), apply(datNP[,11:13], 2, mean), lwd = lwd);
lines(c(0, 1, 2), apply(datP[,11:13], 2, mean), col = 'darkgrey', lwd = lwd);
lines(c(0, 1, 2), apply(datMatched[,10:12], 2, mean), col = 'red', lwd = lwd);
lines(c(0, 1, 2), apply(datNPnomatch[,11:13], 2, mean), col = 'blue', lwd = lwd);
lines(c(0, 1, 2), apply(datNPmatch[,11:13], 2, mean), col = 'green', lwd = lwd);

x = 1.5;
text(x, 2.8, 'NP');
text(x, 2.6, 'P', col = 'darkgrey');
text(x, 2.4, 'Age matched', col = 'red');
text(x, 3, 'NP older', col = 'blue');
text(x, 2.2, 'NP age matched', col = 'green');

# text(x, mean(datNP[, 13]), 'NP');
# text(x, mean(datP[, 13]), 'P', col = 'darkgrey');
# text(x, mean(datMatched[, 12]), 'Age matched', col = 'red');
# text(x, mean(datNPnomatch[, 13]), 'NP older', col = 'blue');
# text(x, mean(datNPmatch[, 13]), 'NP age matched', col = 'green');

