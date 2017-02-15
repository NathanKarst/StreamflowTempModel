files="$(ls *.bil)"

for file in $files
do
	newfile=${file%.*}
	echo $newfile
	gdalwarp -cutline /Users/daviddralle/Dropbox/research/streamflow_temp/model_development/raw_data/projection/bbox.shp -crop_to_cutline -t_srs /Users/daviddralle/Dropbox/research/streamflow_temp/model_development/raw_data/projection/projection.wkt $file $newfile.tif
done

rm *.prj
rm *.csv
rm *.bil
rm *.hdr
rm *.xml 
rm *.stx
rm *.txt
