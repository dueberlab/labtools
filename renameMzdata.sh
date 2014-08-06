for filename in *.mzdata.xml
do
	tr -d '\r' < $filename > tempxml
	rm $filename
	mv tempxml $filename
	newname=$(grep "sampleName" $filename | sed -e 's/[a-zA-Z0-9 <]*\>\([a-zA-Z0-9 _.]*\)\<\/[a-zA-Z0-9 >]*/\1.mzdata.xml/' -e 's/ /_/' -e 's/ /_/' -e 's/ /_/' | sed -e 's/\r//g')
	mv $filename $newname
done
