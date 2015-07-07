setwd('~/Documents/_Fernald_lab');

db = read.csv('Database 05-02-13.csv', na.strings = '');
names(db);
names(db)[c(4, 5, 15)] = c('WEIGHT_FISH', 'LENGTH_FISH', 'WEIGHT_TISSUE');
names(db) = toupper(names(db));

fishIDs = names(table(db$FISH_ID));
dblist = list();
for (id in 1:length(fishIDs))
{
	dblist[[id]] = db[db$FISH_ID == fishIDs[id], ];
	names(dblist)[id] = fishIDs[id];
}
rm(id);

fishWeights = c();
for (id in 1:length(dblist))
{
	w = unique(dblist[[id]]$WEIGHT_FISH);
	if (length(w) == 1)
	{
		fishWeights = c(fishWeights, w);
		names(fishWeights)[id] = names(dblist)[id];
	}
	else 
	{
		cat('fish ', names(dblist)[id], ' has more than one weight, taking mean\n', sep = '');
		fishWeights = c(fishWeights, mean(dblist[[id]]$WEIGHT_FISH));
		names(fishWeights)[id] = names(dblist)[id];
	}
}
rm(id, w);