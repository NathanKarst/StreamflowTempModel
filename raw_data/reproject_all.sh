files="$(ls *.bil)"

for file in $files
do
	newfile=${file%.*}
	echo $newfile
	gdalwarp -cutline ../../projection/bbox.shp -crop_to_cutline -t_srs ../../projection/projection.wkt $file $newfile.tif
done

rm *.prj
rm *.csv
rm *.bil
rm *.hdr
rm *.xml 
rm *.stx
rm *.txt
