cv <-
function(data)
{
	value = sd(as.numeric(data), na.rm = T) / mean(as.numeric(data), na.rm = T);
	return(value);	
}
