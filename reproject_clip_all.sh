DIR=$1

directories="$(ls -d $DIR*)"

for directory in $directories
do
	current_directory="$(pwd)"
	cd $directory
	files="$(ls *.tif)"
	for file in $files
	do
		newfile=${file%.*}
		echo $newfile
		gdalwarp -cutline ../../watershed_poly/eel.shp -crop_to_cutline $file ./blah.tif
		rm $file
		mv ./blah.tif $newfile.tif
	done
	cd $current_directory
done




