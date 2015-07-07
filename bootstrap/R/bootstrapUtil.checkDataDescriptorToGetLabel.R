bootstrapUtil.checkDataDescriptorToGetLabel <-
function (dataDescriptor)
{
	if (!is.null(dataDescriptor) & is.character(dataDescriptor)) {data.lab = dataDescriptor}
	else {
		data.lab = '';
		warning('Please label your data using arg \'dataDescriptor\'\n  It\'s better for everyone');
	}
	return(data.lab);
}
