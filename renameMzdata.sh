for filename in *.mzdata.xml
do
	newname=$(grep "sampleName" $filename | sed -e 's/[a-zA-Z0-9 <]*\>\([a-zA-Z0-9 _.]*\)\<\/[a-zA-Z0-9 >]*/\1.mzdata.xml/' -e 's/ /_/' -e 's/ /_/' -e 's/ /_/' | tr -d '\r' | sed -e 's/\r//g')
	mv $filename $newname
done
