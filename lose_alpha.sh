
files="$(ls ./*.tif)"

for file in $files
do
	newfileappend='temp'
	newfile=${file%.*}
	newfile=$newfile$newfileappend
	gdal_translate -b 1 $file $newfile.tif
	rm $file
	tiffappend=$'.tif'
	newfile=$newfile$tiffappend
	mv $newfile $file
done
	

